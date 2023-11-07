function [FromCam2W, XYZ, RGB] = rigid_transforms(imgseq,w_frame,cam_params, max_n_points)

    FromCam2W(w_frame).R = eye(3);
    FromCam2W(w_frame).T = zeros(3,1);

    % Numero de imagens
    n = size(imgseq,2);
    
    for m = 1:n-1
        
        %Nome das imagens
        depth1_name = imgseq(m).depth;
        depth2_name = imgseq(m+1).depth;
        rgb1_name = imgseq(m).rgb;
        rgb2_name = imgseq(m+1).rgb;

        %Valores RGB para saber no final
        rgb1_values = imread(rgb1_name);
        rgb2_values = imread(rgb2_name);

        %Valores em gray scale para as imagens 1 e 2
        rgb1 = rgb2gray(im2single(rgb1_values));
        rgb2 = rgb2gray(im2single(rgb2_values));

        %Valores de depth       
        depth1 = imread(depth1_name);
        depth2 = imread(depth2_name);
        
        %Matrizes das câmaras
        Kdepth = cam_params.Kdepth;
        Krgb = cam_params.Krgb;

        %Matrizes das câmaras        
        Rdrgb=cam_params.R;
        Tdrgb=cam_params.T;

        %valores xyz e rgb em coordenadas 3D no referencial da camara  
        xyz1_full = get_xyzasus(depth1(:), size(rgb1), 1:size(rgb1,1)*size(rgb1,2), Kdepth,1,0);
        xyz2_full = get_xyzasus(depth2(:), size(rgb2), 1:size(rgb2,1)*size(rgb2,2), Kdepth,1,0);
        xyz1_full(isnan(xyz1_full)) = 0;
        xyz2_full(isnan(xyz2_full)) = 0;
        rgbd1 = get_rgbd(xyz1_full, rgb1_values, Rdrgb, Tdrgb, Krgb);
        rgbd2 = get_rgbd(xyz2_full, rgb2_values, Rdrgb, Tdrgb, Krgb);

        %SIFT
        [coord1, d1] = vl_sift(rgb1);
        [coord2, d2] = vl_sift(rgb2);
        [matches, ~] = vl_ubcmatch(d1, d2) ;

        %coordenadas em pixeis das matches 
        x1 = round(coord1(1:2,matches(1,:))); 
        x2 = round(coord2(1:2,matches(2,:)));

        %Remover matches com depth = 0 e depth > 4000 porque são pontos sem informação
        [xyz1,xyz2,x1] = filter_matches(x1,x2,depth1,depth2,Kdepth);

        %RANSAC
        bestInliersList = ransac(x1,xyz1,xyz2,n);
        
        %Calcular R,T finais com os inliers determinados com o RANSAC
        [R,T] = get_transformation(xyz1,xyz2,bestInliersList);
        
        if m < w_frame        
            %Transformacoes do frame f para min([w_frame m])  
            %inv(R)=R' (R ortogonal)
            for f = 1:1:min([w_frame m])          
                if f==m 
                    FromCam2W(f).R = R';
                    FromCam2W(f).T = -R'*T;    
                else
                    FromCam2W(f).R = R'*FromCam2W(f).R;
                    FromCam2W(f).T = R'*FromCam2W(f).T-R'*T;
                end
            end
        elseif m >= w_frame
            %Transformacoes do frame m+1 para o w_frame 
            FromCam2W(m+1).R=FromCam2W(m).R*R;
            FromCam2W(m+1).T=FromCam2W(m).R*T+FromCam2W(m).T;
        end

        %reshape para um vector
        rows = size(rgb1,1)*size(rgb1,2);
        rgbd1 = reshape(rgbd1, [rows, 3]);
        rgbd2 = reshape(rgbd2, [rows, 3]);

        %inicializacao das vaiaveis dos pontos
        if m==1
            rgb_final=rgbd1;
            xyz_camera{1}=xyz1_full';
        end

        %guardar coordenadas dos pontos no referecial da camara m+1
        rgb_final=[rgb_final; rgbd2];
        xyz_camera{m+1}=xyz2_full';
    end
    
    %converter os pontos do referecial da camara m para o w_frame
    xyz_final=[];
    for m = 1:n
        
        RT=[FromCam2W(m).R,FromCam2W(m).T];
        
        %Aplicar transformação aos pontos da imagem m (em coordenadas homogeneas)
        xyz = RT*[xyz_camera{m}; ones(1,size(xyz_camera{m},2))];
        xyz_final=[xyz_final xyz];
    end
   
    %mostrar point cloud e escalar a pc para max_n_points
    pc = pointCloud(xyz_final','Color',rgb_final);

    if max_n_points < pc.Count
        percentage = max_n_points/pc.Count;
        pc = pcdownsample(pc, 'random' , percentage);
    end
     
    %variaveis de saida
    XYZ=pc.Location;
    RGB=pc.Color;
end


function [xyz1,xyz2,x1] = filter_matches(x1,x2,depth1,depth2,Kdepth)

    i = 1;
    nmatches = size(x1,2);
    while i < nmatches
        positive_depth = depth1(x1(2,i),x1(1,i)) ~= 0 && depth2(x2(2,i),x2(1,i)) ~= 0;
        valid_dist = depth1(x1(2,i),x1(1,i)) < 4000 && depth2(x2(2,i),x2(1,i)) < 4000;
        if  positive_depth && valid_dist 

            x1(3,i) = double(depth1(x1(2,i),x1(1,i)))/1000;
            x1(1:2,i) = x1(1:2,i)*x1(3,i);
            x2(3,i) = double(depth2(x2(2,i),x2(1,i)))/1000;
            x2(1:2,i) = x2(1:2,i)*x2(3,i);
            i = i + 1;
        else
            x1(:,i) = [];
            x2(:,i) = [];
            nmatches = nmatches - 1;
        end
    end
    
    %coordenadas 3d das matches
    xyz1 = Kdepth\x1;
    xyz2 = Kdepth\x2;
end
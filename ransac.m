function [bestInliersList] = ransac(x1,xyz1,xyz2,n)

    % Aproximar numero de iteracoes  
    if n <= 5
        w = 0.4;
    elseif n > 5 && n <= 10
        w = 0.45;
    elseif n > 10
        w = 0.6;
    end
    n_iter = round(log10(1-0.999)/log10(1-w^4));
    
    %RANSAC
    bestInliers = 0;
    for i = 1: n_iter

        indexes = randi(size(x1,2),1,4);

        random_points1 = ones(3,4);
        random_points2 = ones(3,4);
        a = 1;
        for j = indexes %Escolher pontos aleatórios
            random_points1(:,a) = xyz1(:,j);
            random_points2(:,a) = xyz2(:,j);
            a = a + 1;
        end

        %Estimativa do modelo RT desta iteração do RANSAC
        [random_points1, centroid1] = calculateAndSubtractCentroid(random_points1);
        [random_points2, centroid2] = calculateAndSubtractCentroid(random_points2);
        [U, ~, V] = svd(random_points1*random_points2');
        R = U*V';
        T = centroid1 - R*centroid2;
        xyz2_transformed = R*xyz2 + T*ones(1,size(xyz2, 2));

        %Calcular distâncias entre as matches transformados e as matches
        %da outra imagem para saber quais são inliers. O threshold 0.2
        %representa os 20cm de tolerância do Kinect porque as coordenadas
        %xyz estão em metros.
        dists = ones(1,size(xyz1,2)); %Alocar
        for a = 1:size(xyz1,2)
            x_dif =  xyz1(1,a) - xyz2_transformed(1,a); %Distância x
            y_dif =  xyz1(2,a) - xyz2_transformed(2,a); %Distância y
            z_dif =  xyz1(3,a) - xyz2_transformed(3,a); %Distância z
            dists(a) = sqrt( x_dif.^2 + y_dif.^2 +z_dif.^2 );%Distância 3d
        end

        inliers = sum(dists < 0.2); %Índices das matches com erro menor que 20cm

        if (inliers >= bestInliers)

            bestInliersList = (dists < 0.2); %Guardar que pontos são inliers
            bestInliers = inliers; %Guardar o melhor número de inliers obtido
        end   
        
        bestInliersList(2,:)=bestInliersList(1,:);
        bestInliersList(3,:)=bestInliersList(1,:);
    end
end
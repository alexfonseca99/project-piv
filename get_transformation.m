% Calculo de RT atraves da rresolucao do problema de procrustes
function [R,T] = get_transformation(xyz1,xyz2,bestInliersList)

xyz1 = bestInliersList.*xyz1; %Se não for inlier as coordenadas ficam 0
xyz1(:,all(xyz1 == 0))=[]; %Tirar colunas com coordenadas a 0
xyz2 = bestInliersList.*xyz2;
xyz2(:,all(xyz2 == 0))=[];

%Obter centroide e as coordenadas normalizadas
[xyz1,centroid1] = calculateAndSubtractCentroid(xyz1);
[xyz2,centroid2] = calculateAndSubtractCentroid(xyz2);

%Calcular R e T
[U,~,V] = svd(xyz1*xyz2');
R = U*V';
T = centroid1 - R*centroid2;
end
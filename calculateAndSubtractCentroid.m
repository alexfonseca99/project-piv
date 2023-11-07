%pointsIn deve ser 3xn
%Calcula e subtrai o centroide dos pontos pointsIn
function [pointsOut, centroid] = calculateAndSubtractCentroid(pointsIn)
    n = size(pointsIn, 2);
    centroid_x = sum(pointsIn(1,:))/n;
    centroid_y = sum(pointsIn(2,:))/n;
    centroid_z = sum(pointsIn(3,:))/n;
    centroid = [centroid_x; centroid_y; centroid_z];
    pointsOut = pointsIn - centroid*ones(1,n);
end
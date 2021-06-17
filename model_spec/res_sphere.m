function [ d ] = res_sphere( X,M )
%RES_PLANE
model = M(:);
points = X';
d = abs(sqrt((points(:,1)-model(1)).^2 + (points(:,2)-model(2)).^2 + ...
    (points(:,3)-model(3)).^2) - model(4));
d = d(:);
end


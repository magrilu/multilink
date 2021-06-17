function par = fit_sphere(X)
num_points = size(X,2);
if(num_points==4)
    par =  fit_sphere_cardmss(X);
elseif(num_points >4)
    [c,r] = sphereFit(X');
    par = [c(1);c(2);c(3);r];
end
par = par(:);
end




function model = fit_sphere_cardmss(points)
points = points(1:3,:)';
X = [points, ones(size(points,1),1)];
m11 = det(X);
if abs(m11)<=eps(class(points))
    model = cast([], 'like', points);
    return;
end

X(:,1) = points(:,1).^2+points(:,2).^2+points(:,3).^2;
m12 = det(X);

X(:,2) = X(:,1);
X(:,1) = points(:,1);
m13 = det(X);

X(:,3) = X(:,2);
X(:,2) = points(:,2);
m14 = det(X);

X(:,1) = X(:,3);
X(:,2:4) = points;
m15 = det(X);

a =  0.5*m12/m11;
b =  0.5*m13/m11;
c =  0.5*m14/m11;
d = sqrt(a^2+b^2+c^2-m15/m11);
model = cast([a, b, c, d], 'like', points);
end
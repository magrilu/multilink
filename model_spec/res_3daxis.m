function d = res_3daxis(X,M)
%RES_3DLINE 3D axis is defined as a pair of point on it M = <a,b>
% and the distance from <a,b> is calculated as 
% http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
n = size(X,2);
x1 = M(1:3);
x2 = M(4:6);
v = x2-x1;
den = norm(v);
w = x1-X;
d = sqrt(sum(cross(repmat(v,1,n),w).^2,1))./den;
end


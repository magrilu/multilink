function model = fit_pc_cylinder(X)

p1 = X(1:3,1)';
n1 = X(4:6,1)';
p2 = X(1:3,2)';
n2 = X(4:6,2)';
w = p1 + n1 - p2;


a = n1 * n1';
b = n1 * n2';
c = n2 * n2';
d = n1 * w';
e = n2 * w';
D = a * c - b * b;
if abs(D) < 1e-5
    s = 0;
    if b > c
        t = d / b;
    else
        t = e / c;
    end
else
    s = (b * e - c * d) / D;
    t = (a * e - b * d) / D;
end

% P0 is a point on the axis
p0 = p1 + n1 + s * n1;
% dp is a normalized vector for the direction of the axis
dp = p2 + t * n2 - p0;
dp = dp / norm(dp);

% figure; hold all;
% plot3(p1(1),p1(2),p1(3));
% plot3(p2(1),p2(2),p2(3));
% quiver3(p1(1),p1(2),p1(3),n1(1),n1(2),n1(3),'k')
% quiver3(p2(1),p2(2),p2(3),n2(1),n2(2),n2(3),'k')
% plot3(p0(1),p0(2),p0(3));
% quiver3(p0(1),p0(2),p0(3),dp(1),dp(2),dp(3),'c')

p1p0 = p1 - p0;
p2p1 = dp;
c = [p1p0(2)*p2p1(3) - p1p0(3)*p2p1(2), ...
     p1p0(3)*p2p1(1) - p1p0(1)*p2p1(3), ...
     p1p0(1)*p2p1(2) - p1p0(2)*p2p1(1)];
% p2p1 is a unit vector, so the denominator is not needed 
r = sqrt(sum(c.*c, 2));

model = [p0, dp, r];
model = model';

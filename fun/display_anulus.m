function [  ] = display_anulus(X, par, epsi , col, alphaValue)
%DISPLAY_ANULUS
if(nargin==2)
    epsi=0;
    col = 'k';
    alphaValue = 0.5;
elseif(nargin==3)
     col = 'k';
     alphaValue = 0.5;
elseif(nargin ==4)
    alphaValue = 0.5;
end
radius_out = par(3) +epsi; 
radius_inn = par(3) -epsi;

% % compute the angular range of points
% n = size(X,2);
% U = X - repmat(par(1:2),1,n); % vector from the center to the points
% U = normc(U); % transform in unitary vector
% thetas = nan(n,1);
% for i =1:n
%     thetas(i) = atan2d(U(2,i),U(1,i));
% end
% min_theta = min(thetas);
% max_theta = max(thetas);

if(radius_out>1e3)
    k = 1000; 
else
    k = 100;
end
t = linspace(0, 2*pi,k);
x_inn = radius_out*cos(t) + par(1);
y_inn = radius_out*sin(t) + par(2);
x_out = radius_inn*cos(t) + par(1);
y_out = radius_inn*sin(t) + par(2);
px = [x_inn(1), x_out,flip(x_inn(2:end))];
py = [y_inn(1), y_out,flip(y_inn(2:end))];
patch(px, py, col, 'FaceAlpha',alphaValue,'EdgeColor','None')
%rectangle('Position',[par(1) - radius_inn, par(2) - radius_inn, radius_inn*2, radius_inn*2],'Curvature',[1,1],'FaceColor','w','EdgeColor',col);
end

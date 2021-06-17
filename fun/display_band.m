function [ ] = display_band( X, par,epsi,col,alphaValue )
%DISPLAY_BAND Summary of this function goes here
%   Detailed explanation goes here
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
%%
x1 = min(X(1,:));
x2 = max(X(1,:));
y1 = min(X(2,:));
y2 = max(X(2,:));
% extend the line for the k% of its lenght 
k =0.01;
dx = x2-x1;
dy = y2-y1;
l = sqrt(dx^2+dy^2);
x1 = x1 - k*l;
x2 = x2 + k*l;
y1 = y1 - k*l;
y2 = y2 + k*l;
%%
par = par./norm(par);
if(abs(par(2))<2e-2 )
    %disp('vertical line')
    a = [(par(2)*y1 + par(3))/-par(1); y1];
    b = [(par(2)*y2 + par(3))/-par(1); y2];
else
    a = [x1; (par(1)*x1 + par(3))/-par(2)];
    b = [x2; (par(1)*x2 + par(3))/-par(2)];
end
if(abs(par(1))<1e-5)
    v = [0;1];
else
    v = [1; par(2)/par(1)]; % direction perpendicular to the line;
end  
    v = normc(v);
    % corners of the bands
    a1 = a - epsi.*v;
    a2 = a + epsi.*v;
    b1 = b - epsi.*v;
    b2 = b + epsi.*v;
    
    px = [a1(1),b1(1),b2(1),a2(1)];
    py = [a1(2),b1(2),b2(2),a2(2)];
    patch(px,py,col,'FaceAlpha',alphaValue,'EdgeColor','None');

end


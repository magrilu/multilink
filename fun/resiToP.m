function [P] = resiToP(R,epsi,useX84)
%
if(nargin <3)
    useX84 = false;
end
P = zeros(size(R));
heta = 0.05; % minumum preference for residual equal to epsilon
if(useX84)
    for j= 1:size(R,2)
        inliers = R(:,j)<epsi;
        epsix84 = min(x84_( R(inliers,j), 2),epsi);
        sigma2 = -epsix84^2/log(heta);
        inliers = R(:,j)<=epsix84;
        P(inliers,j) = exp(-(R(inliers,j).^2)./(sigma2));
    end
else
sigma2 = -epsi^2/log(heta);
I = R<epsi;
P(I) = exp(-(R(I).^2)./(sigma2));
end


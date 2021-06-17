function [ mindist ] = distPointParabola(P,H)
%DISTPOINTPARABOLA 
% P point
% H parabola 

gamma = H(3)-P(2); %questa è un'immagine
% derivative
f = [4*H(1)^2, 6*H(1)*H(2), 2*( 1 +H(2)^2 + 2*H(1)*gamma), -2*P(1) + 2*H(2)*gamma ]; % questo diventa un vettore 3D size(immagine)x f f_H(p1,p2)
u = roots(f);


k = numel(u); % number of solutions
Q = nan(2,k); % possible extrema points on the parabola
d = Inf(1,k);

for i=1:k
   if(isreal(u(i)))  
       Q(:,i) = [u(i);H(1)*u(i)^2+H(2)*u(i)+H(3)];
       d(i)   = norm(Q(:,i)-P);
   end
end
[mindist,v] = min(d);

end


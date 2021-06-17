function  d = F_sampson_distance_sqr(F, m1, m2);

% sampson(F, m1, m2)  valuate the first order approximation of the geometric error
% of the fit of F  with respect to a set of matched points m1, m2.
%
% Returns an approximation of the **squared** geometric distance from the
% 4d joint-space point [m1;m2] to the F manifold. 
%
% Reference: Luong, Faugeras. IJCV 1996 (called 'Gradient criterion")  

% Author: A. Fusiello. Extracted and adapted from a piece of code by Peter Kovesi 
% The interface is the same as vgg_H_sampson_distance_sqr (see)

 
[rm1,cm1]=size(m1);
if (rm1 ~= 3)
    error('This function requires homogeneous coordinates');
end

[rm1,cm1]=size(m2);
if (rm1 ~= 3)
    error('This function requires homogeneous coordinates');
end


lng = size(m1,2);

m2tFm1 = zeros(1,lng);
for n = 1:lng
    m2tFm1(n) = m2(:,n)'*F*m1(:,n);
end

Fm1 = F*m1;
Ftm2 = F'*m2;

% Evaluate quared distances
d =  m2tFm1.^2 ./ (Fm1(1,:).^2 + Fm1(2,:).^2 + Ftm2(1,:).^2 + Ftm2(2,:).^2);



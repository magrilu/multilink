function [par] = fit_line_from_circle( X, center )
% define a line perpendicular

% angular coefficient of the line passing through X and the center
den = (X(1)-center(1));
num = (X(2)-center(2));
m1 =  num/den;
if(abs(den)<1e-4)
    % line is vertical
    m2 = 0;
else
    m2 = -1/m1;
end
   par = fit_line_parallel(X,m2);
end
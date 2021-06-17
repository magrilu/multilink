function [par] = fit_line_parallel( X, m )
if(isinf(m))
    %line vertical
    par = [1; 0; -X(1)];
else
    q = X(2) - m*X(1);
    par = [m; -1; q];
end
par = par./(norm(par(1:2)));
end
function d  = res_line( X, M )
%RES_LINE compute the residual between points and a line
num_lines = size(M,2);

if(num_lines==1)
    d =  abs(M(1)*X(1,:) + M(2)*X(2,:) + M(3));
    d = d(:);
else
    n = size(X,2);
    d = nan(n, num_lines);
    for i = 1:num_lines
         d(:,i) =  abs(M(1,i)*X(1,:) + M(2,i)*X(2,:) + M(3,i));
    end
end


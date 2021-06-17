function  d  = res_circle( X, M )
%RES_CIRCLE compute the residual between points and circle
n = size(X,2);
num_circles = size(M,2);
if(num_circles==1)
 d = abs( sqrt(sum((X-repmat(M(1:2),1,n)).^2,1))-M(3) );
 d = d(:);
else
    d = nan(n,num_circles);
    for i = 1:num_circles
        d(:,i) = abs( sqrt(sum((X-repmat(M(1:2,i),1,n)).^2,1))-M(3,i) );
    end

end
end


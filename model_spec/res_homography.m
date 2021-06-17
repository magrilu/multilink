function  d  = res_homography( X, M )
%RES_CIRCLE compute the residual between point correspondences and
%fundamental matrix
n = size(X,2);
num_fm = size(M,2);
m1 = X(1:3,:);
m2 = X(4:6,:);

d = nan(n,num_fm);
for i = 1:num_fm
    H = reshape(M(:,i),[3,3]);
  
    if(0)
        % geometric
        hm1 = homo2cart(H*m1);
        m2h = homo2cart(H\m2);
        e1 = sum((m1 - m2h).^2);
        e2 = sum((m2 - hm1).^2);
        d(:,i) = e1+ e2;
    else
        % sampson
        d(:,i) = sqrt(vgg_H_sampson_distance_sqr(H,m1,m2));
    end
    
    
end

end


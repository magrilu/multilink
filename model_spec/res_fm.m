function  d  = res_fm( X, M )
%RES_CIRCLE compute the residual between point correspondences and
%fundamental matrix
n = size(X,2);
num_fm = size(M,2);
m1 = X(1:3,:);
m2 = X(4:6,:);

    d = nan(n,num_fm);
    for i = 1:num_fm
        F = reshape(M(:,i),[3,3]);      
%         Fm1 = F*m1;
%         l2 = Fm1./ sqrt(sum(Fm1(1:2,:).^2)); % normalized epiploar lines in the second image
%         Fm2 = F'*m2;
%         l1 = Fm2./ sqrt(sum(Fm2(1:2,:).^2)); % normalized epiploar lines in the first image
%         res1 = sqrt(sum(abs(l1.*m1)));
%         res2 = sqrt(sum(abs(l2.*m2)));
        d(:,i) = sqrt(F_sampson_distance_sqr(F,m1,m2));
    end

end


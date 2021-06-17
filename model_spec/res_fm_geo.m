function  d  = res_fm_geo( X, M )

n = size(X,2);
num_fm = size(M,2);
m1 = X(1:3,:);
m2 = X(4:6,:);

    d = nan(n,num_fm);
    for i = 1:num_fm
        F = reshape(M(:,i),[3,3]);      
        Fm1 = F*m1;
        l2 = Fm1./ sqrt(sum(Fm1(1:2,:).^2)); % normalized epiploar lines in the second image
        Fm2 = F'*m2;
        l1 = Fm2./ sqrt(sum(Fm2(1:2,:).^2)); % normalized epiploar lines in the first image
        res1 = abs(sum(l1.*m1));
        res2 = abs(sum(l2.*m2));
       d = res1 + res2;
    end

end


function [flg] = is_homography_from_fund_degen(X, theta,F, inds)
%IS_HOMOGRAPHY_DEGEN 

if(any(isnan(theta))|| any(isinf(theta)))
    flg = true;
    return;
end

x1 = X(1:3,inds);
x2 = X(4:6,inds);

flg1 = iscolinear(x1(:,1),x1(:,2),x1(:,3))|| iscolinear(x2(:,1),x2(:,2),x2(:,3));
if(flg1)
    flg = flg1;
    return;
end
flg2 = validateTheta_homography(X, theta, inds);
H = reshape(theta,[3,3]);
J = H'*F +F'*H;
flg3 = all(J(:)>1e-5);
flg = flg1 || ~flg2 || flg3;
end


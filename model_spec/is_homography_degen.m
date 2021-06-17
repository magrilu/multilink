function [flg] = is_homography_degen(X, theta, inds)
%IS_HOMOGRAPHY_DEGEN 

flg1 = validateMSS_homography(X, inds);
flg2 = validateTheta_homography(X, theta, inds);
flg = ~flg1 || ~flg2;
end


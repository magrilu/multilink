function h = fit_homography(mss)
%FIT_homography 
 %h = HomographyDLT(mss(1:3,:),mss(4:6,:));
% h =  homography2d(mss(1:3,:),mss(4:6,:));
if(size(mss,1)>4)
    mss = remove_repeated_points(mss);
    h = homog_lin(mss(1:2,:),mss(4:5,:));
else
    h = homog_lin(mss(1:2,:),mss(4:5,:));
end
 h = h(:);
end


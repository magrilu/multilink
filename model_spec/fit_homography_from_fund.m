function h = fit_homography_from_fund(mss,A,e2)
%FIT_homography 
x1 = mss(1:3,:);
x2 = mss(4:6,:);
M = [x1'];
if(rcond(M)<1e-10)
    h = nan(9,1);
    return;
end
b = nan(3,1);
for i = 1:3
    den = cross(x2(:,i),e2);
    b(i) = (cross(x2(:,i),A*x1(:,i))'* den)./norm(den).^2;
end
H = A - e2*(M\b)'; 
h = H(:);

end


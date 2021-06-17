function flg =  is_cylinder_degen(mss)
% flg = true if the cylinder is degenrate
tol = 1e-4;
%chek on the position of points
x1 = mss(1:3,1);
x2 = mss(1:3,2);
if(norm(x1-x2)< tol )
    flg1 = true;
    mss
else
    flg1 = false;
end
% check on the normal
n1 = mss(4:6,1);
n2 = mss(4:6,2);
s = cross(n1,n2);
if(norm(s)<tol)
     flg2 = true;
     mss
else
     flg2 = false;
end

flg = flg1 || flg2;
end
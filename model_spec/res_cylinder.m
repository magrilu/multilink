function [ d ] = res_cylinder( X, M )
%RES_CYLINDER
% M represent a cylinder M(1:3) direction W must be normalized
%                        M(4:6) point on the axis
%                        M(7)   radius

num_cylinder = size(M,2);

if (num_cylinder ==1)
    W  = M(1:3);
    PC = M(4:6);
    r = M(7);
    P = eye(3) - W*W';
    XC = bsxfun(@minus,X(1:3,:),PC);
    res_geo = abs(sqrt(sum( (XC'*P)'.* XC))-r);
    d= res_geo;
else
    for i = 1:num_cylinder
        W  = M(1:3,i);
        PC = M(4:6,i);
        r = M(7,i);
        P = eye(3) - W*W';
        XC = bsxfun(@minus,X(1:3,:),PC);
        res_geo = abs(sqrt(sum( (XC'*P)'.* XC))-r);
        d(:,i)= res_geo;
        
    end
end


end


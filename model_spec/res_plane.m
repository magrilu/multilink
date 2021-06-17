function [ d ] = res_plane( X,M )
%RES_PLANE
num_planes = size(M,2);
if(num_planes==1)
    d = abs(M(1)*X(1,:) + M(2)*X(2,:) + M(3)*X(3,:) + M(4));
else
    n = size(X,2);
    d = nan(n, num_planes);
    for i = 1:num_planes
        d(:,i) =  abs(M(1,i)*X(1,:) + M(2,i)*X(2,:) + M(3,i)*X(3,:)+ M(4,i));
    end
end


end


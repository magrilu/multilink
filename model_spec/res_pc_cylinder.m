function  d = res_pc_cylinder( X, M )
%RES_CYLINDER

n = size(X,2);
m = size(M,2);
d = nan(n,m);

for i = 1:m
    
    % compute the distance between X and the 3Daxis of the cylinder
    % represented by its two supporting points (FINITE MODEL)
    t = res_3daxis(X(1:3,:), M(1:6,i));
    d(:,i) = abs(t -  M(7,i));
    
    %
    %figure;
    %scatter3(X(1,:),X(2,:),X(3,:),10,d(:,i)); hold all;
    %plot3(M(1,i),M(2,i),M(3,i),'sr');
    %plot3(M(4,:),M(5,:),M(6,:),'>g');
    
end
end


function [ Y,flg ] = remove_repeated_points(X)
%REMOVE_REPEATED_POINTS remove data points removing repeated point
%   Detailed explanation goes here
N=size(X,2);
flg=true(1,N);
D=pdist(X','euclidean');
D=squareform(D);
for i=1:N
    for j=i+1:N
        if(D(i,j)==0)
            flg(j)=false;
            %disp('repeated point')
            %X(:,i)
            %X(:,j)
        end
    end
end
Y=X(:,flg);
end


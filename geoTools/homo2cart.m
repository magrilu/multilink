function [ X ] = homo2cart( Q )
%HOMO2CART 
d = size(Q,1);
X= [Q(1:d-1,:)./repmat(Q(end,:),d-1,1); ones(1,size(Q,2))];

end


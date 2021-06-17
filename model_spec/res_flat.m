function resi = res_flat(X,flat)
%RES_FLAT given a set of points X and an affine subspaces
% (flats) compute the point-flat distances and arrange in a matrix d
% Luca Magri
% November 2020 - potrebbe rimpiazzare tutte quelle vecchie

if(isstruct(flat))
    orthProjection = flat.orthProjection;
    origin = flat.origin;
else
    dime = size(X,1);
    origin = flat(1:dime);
    orthProjection = flat(dime+1:end);
    orthProjection = reshape(orthProjection,[dime,dime]);
end

n = size(X,2);
resi = nan(n, 1);
Y = X-origin;
temp = orthProjection * Y;
resi = sqrt(sum(temp.^2,1));
% 
% for i = 1:n
% resi(i) =  norm(orthProjection * Y(:,i));
% end


end


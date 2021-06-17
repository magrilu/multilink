function [flat] = fit_flat(X,dimflat)
% FIT_FLAT
% Fit an affine subspace (flat) of given dimension (dimflat) to a set of
% points X.
% INPUT:
% X d-dimensional points organized by columns
% dimFlat the dimension of the flat (it should be less than d)
% OUTPUT:
% flat is a structure whose fields are:
% - origin is a point on the affine subspace
% - basis is a basis for the direction subspace of the affine subspace (la giacitura)
% - orthProjection is matrix that return the orthogonal components of a point -
%   origin with respect to the orthogonal subspace of the basis:
%   this matrix can be used to compute the distance of a point P to the flat:
%   orthProjection * (P-origin) is the vector from P to the projection of P
%   onto the flat, see res_flat.
%
% References:
% - https://www.geometrictools.com/Documentation/LeastSquaresFitting.pdf
% - https://ocw.mit.edu/courses/mathematics/18-06sc-linear-algebra-fall-2011/least-squares-determinants-and-eigenvalues/projections-onto-subspaces/MIT18_06SCF11_Ses2.2sum.pdf
%
% Luca Magri
% November 2020 - potrebbe rimpiazzare tutte quelle vecchie

assert(size(X,2)>=dimflat+1,'Not enough points to fit an affine subspace');

origin = mean(X,2);
Y = X-origin; % centered data

[U,~,~] = svd(Y);
B = U(:,1:dimflat);

flat.origin = origin;
flat.basis = B;
flat.orthProjection = eye(size(X,1))-B*B';
end


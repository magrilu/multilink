function [gScore, dataFidelity, modelComplexity] = getGricScore(rSqr,sigma,r,d,k,lambda1, lambda2)
%GETGRICSCORE 
% k number of parameters;
% d dimension of the manifold
% r dimension of the space

n = numel(rSqr);
dataFidelity = sum(min(rSqr./sigma^2, 2*(r-d)));
modelComplexity =  lambda1*n*d + lambda2*k;
gScore = dataFidelity + modelComplexity;



end


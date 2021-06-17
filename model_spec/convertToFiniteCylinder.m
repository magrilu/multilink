function modelParams = convertToFiniteCylinder(modelParams, X)
% Get the direction of cylinder axis
dp = modelParams(4:6);
% Get the inlier points
points = X(1:3,:);
n = size(X,2);
% Get a point on the axis
p0 = modelParams(1:3);
% Describe the axis as a line equation: p0 + k * dp, and find the
% projections of inlier points on this line
%k = sum(points.*repmat(dp, n, 1),2) - p0 * dp';
k = sum((points - p0).* repmat(dp,1,n),1);
% Find the two extreme points
pa = p0 + min(k) * dp;
pb = p0 + max(k) * dp;
% Set to the parameters required by cylinderModel object
modelParams(1:6) = [pa, pb];

% figure; hold all;
% scatter3(X(1,:),X(2,:),X(3,:));
% plot3(pa(1),pa(2),pa(3),'r*');
% plot3(pb(1),pb(2),pb(3),'g*');
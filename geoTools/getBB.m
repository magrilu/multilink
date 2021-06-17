function [bb] = getBB(X,k)
%GETBB get the bounding box of points X
% the bounding box is expanded a little bit, the bigger  k the more it
% is extended
%
%
% to do: multidimensionale


xmin = min(X(1,:));
xmax = max(X(1,:));
ymin = min(X(2,:));
ymax = max(X(2,:));
wx = abs(xmax - xmin);
wy = abs(ymax -ymin);
dx = k * wx;
dy = k * wy;

bb.xmin = xmin - dx;
bb.xmax = xmax + dx;
bb.ymin = ymin - dy;
bb.ymax = ymax + dy;
end


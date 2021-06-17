function [ W, PC, r ,e] = fit_cylinder( X,w0 )
%FIT_CYLINDER 
%
%
% author: Luca Magri 2018
% reference: Fitting 3D data with a Cylinder - David Eberly
%
% C point on axis
% r radius of cylinder
% W axis unit direction
% 6 parameters


% subtract mean
t = mean(X,2);
X = bsxfun(@minus, X,t);


% if the data can be described by a line, you can use the eigenvector of the
% covariance matrix
%[u,v]=eig(X*X');
%[~,t] = max(diag(v));
%W_start = normc(u(:,t));

% use w0 as the starting direction for simulated annealing
if(nargin<2 || isempty(w0))
    s0 = [0,pi/2];
    options = optimoptions(@simulannealbnd,'Display','off');
    
else
    [el,az,~] = cart2sph(w0(1),w0(2),w0(3));
    s0 = [el,az];
    options = optimoptions(@simulannealbnd,'Display','off','MaxIterations',5);
    
end

fun = @(s) cost(s,X);

lb = zeros(2,1);
ub = [2*pi,pi/2];
s = simulannealbnd(fun,s0,lb,ub,options);
W = nan(3,1);
[W(1),W(2),W(3)] = sph2cart(s(1),s(2),1);



[e, PC, r] = G_cylinder(W,X);
 
PC = PC+t;
  
 
   

end

%% 
function [error] = cost(s,X)
% auxiliary function for simulated annealing
[Wx,Wy,Wz] = sph2cart(s(1),s(2),1);
[error, ~, ~] = G_cylinder([Wx;Wy;Wz],X);
end

%%
function [error, PC, r] = G_cylinder(W,X)
%function [error, PC, r] = G_cylinder(W,X)
% evaluates the function G(W) and generate the corresponding PC and r

n = size(X,2);
W = normc(W);
P = eye(3) - W*W';
S = [0 -W(3), W(2); W(3),0, -W(1); -W(2), W(1), 0];
Y = P*X;
A = (Y*Y')./n;
A_hat = S*A*S';
L = sum(Y.*Y,1); % squared Length
aL = mean(L);    % averaged Length
B = mean(bsxfun(@times, L, Y),2);
PC = (A_hat*B)./sum(diag(A_hat*A)); % point on the axis of the cylinder


% radius: mean distance from PC of the projected points
diff = bsxfun(@minus,PC,Y);
norm_diff_sq = sum(diff.*diff,1); % squared distance
r =  sqrt(mean(norm_diff_sq));

% error 

% error = 0;
% for i = 1: n
%     term = L(i)- aL - 2* Y(:,i)'*PC;
%     error = error + term*term; 
% end
% error = error/n;

% algebraic distance 
%res = L-aL -   (2 .*  Y'*PC)';
%error = mean(res.^2);


% geometric distance
 XC = bsxfun(@minus,X,PC);
 res_geo = abs(sqrt(sum( (XC'*P)'.* XC))-r);
 error = mean(res_geo.^2);

end



function [ok, msLC, msOutput] = isMergeableGricPlaneCylinderSphere(X, L, i, j, lambda1 , lambda2, sigma)
% Check if two clusters A and B can be merged.
o_verbose = false;
[okl, msPlane, msOutPlane] = isMergeableGricPlane(X(1:3,:), L, i, j, lambda1 , lambda2, sigma);
[okc, msCylin, msOutCylin] = isMergeableGricCylinder(X, L, i, j, lambda1 , lambda2, sigma);
[oks, msSphere, msOutSphere] = isMergeableGricSphere(X, L, i, j, lambda1 , lambda2, sigma);
%%
scores = [msPlane.gric.before, msPlane.gric.after, msCylin.gric.before, msCylin.gric.after...
    msSphere.gric.before, msSphere.gric.after];
[~, ind] = min(scores);

if(o_verbose)
       f = [msPlane.fidelity.before, msPlane.fidelity.after, msCylin.fidelity.before, msCylin.fidelity.after, msSphere.fidelity.before, msSphere.fidelity.after];
       c = [msPlane.complexity.before, msPlane.complexity.after, msCylin.complexity.before, msCylin.complexity.after, msSphere.complexity.before, msSphere.complexity.after]./lambda2;
       disp('--------------------------------------')
       disp([scores;f;c]);
       disp('--------------------------------------')
end
%%  keep the result with the minimum model selection score   
switch(ind)
    case 1
        % do not merge clusters: seprate lines win
        ok = false;
        msLC = msPlane;
        msOutput = msOutPlane;
    case 2
        % do merge clusters with line
        ok = true;
        msLC = msPlane;
        msOutput = msOutPlane;
        assert(okl);
        
    case 3
        % do not merge clustes: separate circles win
        ok = false;
        msLC = msCylin;
        msOutput = msOutCylin;
    case 4
        % do merge clusters with line
        ok = true;
        msLC = msCylin;
        msOutput = msOutCylin;
        assert(okc);
        
    case 5
        % do not merge clustes: separate circles win
        ok = false;
        msLC = msSphere;
        msOutput = msOutSphere;
    case 6
        % do merge clusters with line
        ok = true;
        msLC = msSphere;
        msOutput = msOutSphere;
        assert(oks);
end




end

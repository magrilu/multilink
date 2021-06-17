function [ok, msLCP, msOutput] = isMergeableGricLCP(X, L, i, j, lambda1 , lambda2, sigma)
% Check if two clusters A and B can be merged.
o_verbose = false;
[okl, msLine, msOutLine] = isMergeableGricLine(X, L, i, j, lambda1 , lambda2, sigma);
[okc, msCirc, msOutCirc] = isMergeableGricCircle(X, L, i, j, lambda1 , lambda2, sigma);
[okp, msPara, msOutPara] = isMergeableGricParabola(X,L,i,j,lambda1, lambda2, sigma);
%%
scores = [msLine.gric.before, msLine.gric.after,...
    msCirc.gric.before, msCirc.gric.after,...
    msPara.gric.before, msPara.gric.after];
[~, ind] = min(scores);

if(o_verbose)
    f = [msLine.fidelity.before, msLine.fidelity.after,...
        msCirc.fidelity.before, msCirc.fidelity.after,...
        msPara.fidelity.before, msPara.fidelity.after];
    c = [msLine.complexity.before, msLine.complexity.after,...
        msCirc.complexity.before, msCirc.complexity.after...
        msPara.complexity.before, msPara.complexity.after]./lambda2;
    disp('--------------------------------------')
    disp([scores;f;c]);
    disp('--------------------------------------')
end
%%  keep the result with the minimum model selection score
switch(ind)
    case 1
        % do not merge clusters: seprate lines win
        ok = false;
        msLCP = msLine;
        msOutput = msOutLine;
    case 2
        % do merge clusters with line
        ok = true;
        msLCP = msLine;
        msOutput = msOutLine;
        assert(okl);
        
    case 3
        % do not merge clustes: separate circles win
        ok = false;
        msLCP = msCirc;
        msOutput = msOutCirc;
    case 4
        % do merge clusters with line
        ok = true;
        msLCP = msCirc;
        msOutput = msOutCirc;
        assert(okc);
        
    case 5
        % do not merge clustes: separate parabolas win
        ok = false;
        msLCP = msPara;
        msOutput = msOutPara;
    case 6
        % do merge clusters with parabola
        ok = true;
        msLCP = msPara;
        msOutput = msOutPara;
        assert(okp);   
end




end

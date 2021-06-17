function [ok, msLCP, msOutput] = isMergeableGricHAF(X, L, i, j, lambda1 , lambda2, sigma)
% Check if two clusters A and B can be merged.
o_verbose = false;
[okl, msHomo, msOutLine] = isMergeableGricHomography(X, L, i, j, lambda1 , lambda2, sigma);
[okc, msAfun, msOutCirc] = isMergeableGricAffineFundamental(X, L, i, j, lambda1 , lambda2, sigma);
[okp, msFund, msOutPara] = isMergeableGricFundamental(X,L,i,j,lambda1, lambda2, sigma);
%%
scores = [msHomo.gric.before, msHomo.gric.after,...
    msAfun.gric.before, msAfun.gric.after,...
    msFund.gric.before, msFund.gric.after];
[~, ind] = min(scores);

if(o_verbose)
    f = [msHomo.fidelity.before, msHomo.fidelity.after,...
        msAfun.fidelity.before, msAfun.fidelity.after,...
        msFund.fidelity.before, msFund.fidelity.after];
    c = [msHomo.complexity.before, msHomo.complexity.after,...
        msAfun.complexity.before, msAfun.complexity.after...
        msFund.complexity.before, msFund.complexity.after];
    disp('--------------------------------------')
    disp([scores;f;c]);
    disp('--------------------------------------')
end
%%  keep the result with the minimum model selection score
switch(ind)
    case 1
        % do not merge clusters: seprate homographies win
        ok = false;
        msLCP = msHomo;
        msOutput = msOutLine;
    case 2
        % do merge clusters with homo
        ok = true;
        msLCP = msHomo;
        msOutput = msOutLine;
        assert(okl);
        
    case 3
        % do not merge clustes: separate fundamentals affine win
        ok = false;
        msLCP = msAfun;
        msOutput = msOutCirc;
    case 4
        % do merge clusters with fundamental affine
        ok = true;
        msLCP = msAfun;
        msOutput = msOutCirc;
        assert(okc);
        
    case 5
        % do not merge clustes: separate fundamentals win
        ok = false;
        msLCP = msFund;
        msOutput = msOutPara;
    case 6
        % do merge clusters with fundamental
        ok = true;
        msLCP = msFund;
        msOutput = msOutPara;
        assert(okp);   
end




end

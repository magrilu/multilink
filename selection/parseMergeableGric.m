function [isMergeableGric,cardmss] = parseMergeableGric(modelStr)
if(strcmp(modelStr,'line'))
    isMergeableGric = @isMergeableGricLine;
    cardmss =2;
elseif(strcmp(modelStr,'circle'))
    isMergeableGric = @isMergeableGricCircle;
    cardmss = 3;
elseif(strcmp(modelStr,'lc'))
    isMergeableGric = @isMergeableGricLC;
    cardmss = 3;
elseif(strcmp(modelStr,'lcp'))
    isMergeableGric = @isMergeableGricLCP;
    cardmss = 3;
elseif(strcmp(modelStr,'homography'))
    isMergeableGric = @isMergeableGricHomography;
    cardmss = 4;
elseif(strcmp(modelStr,'affine_fundamental'))
    isMergeableGric = @isMergeableGricAffineFundamental;
    cardmss = 4;
elseif(strcmp(modelStr,'fundamental'))
    isMergeableGric = @isMergeableGricFundamental;
    cardmss = 8;
elseif(strcmp(modelStr,'haf'))
    isMergeableGric = @isMergeableGricHAF;
    cardmss = 8;
elseif(strcmp(modelStr,'planecylinder'))
    isMergeableGric = @isMergeableGricPlaneCylinder;
    cardmss = 3;
elseif(strcmp(modelStr,'plane'))
    cardmss = 3;
    isMergeableGric = @isMergeableGricPlane;
elseif(strcmp(modelStr,'sphere'))
    isMergeableGric = @isMergeableGricSphere;
    cardmss = 3;
elseif(strcmp(modelStr,'planecylindersphere'))
    isMergeableGric = @isMergeableGricPlaneCylinderSphere;
    cardmss = 4;
end
end
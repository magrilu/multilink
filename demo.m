% This snippet provides a simple demo of Multi-Link
% a multi-class multi-model fitting clustering algorithm based on
% preferences, presented in:  
% % Magri, Leveni, Boracchi, MultiLink: Multi-class Structure Recovery 
% via Agglomerative Clustering and Model Selection. CVPR 2021

%%
addpath(genpath('../multiLink_beta'))
%%  load data
clear all;
close all;
temp = load('data/ticktoe.mat');
X = temp.X;
G = temp.G;
clear temp;

figure;
gscatter(X(1,:),X(2,:),G);
legend off, axis off, axis equal;
title('Input Data');

%% preference embedding with multiple model class

% mixed hypotheses sampling: lines and circles
optsSampling.m = 2000; % number of hypotheses
optsSampling.sampling = 'localized';
optsSampling.robust = 'on';
optsSampling.voting = 'gauss';
% sampling lines
optsl = optsSampling;
optsl.model = 'line';
Sl = computeResi(X,optsl);
% sampling circles
optsc = optsSampling;
optsc.model = 'circle';
Sc =  computeResi(X,optsc);
% preference computation
epsi = 2e-2; % inlier threhsold
[Sl.P] = resiToP(Sl.R,epsi);
[Sc.P] = resiToP(Sc.R,epsi);
P =[Sl.P,Sc.P];

%% agglomerative clustering with model selection
modelType = 'lc'; % alternative  models line (l) and circle (c)
gricParam.lambda1 = 1;
gricParam.lambda2 = 2;
gricParam.sigma = epsi;
C = multiLink(X,P,modelType,gricParam);
thCard = 10; % prune small clusters, e.g. using the cardinality of the mss
Cpruned = prune_small_clust(C,thCard);

%% estimated clusters
figure; 
gscatter(X(1,:),X(2,:),Cpruned);
legend off, axis off, axis equal;
title('Estimated Clusters');
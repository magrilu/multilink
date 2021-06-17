function [] = displaySLMerge(X,L,i,j, n,cmap, img)
%DISPLAYMERGE display the merge during the jslinkage clustering
if(nargin < 7)
    img = [];
end
mrkSize      = 100;
mrkSizeBig   = 3* mrkSize;
mrkSizeClust = 10* mrkSizeBig;
mrkSizePair  = 95;
%scheme = 'Set2';
syms = 'o';'dphv<^>';

%%
set(gcf,'color','w');
if(~isempty(img))
    if(numel(size(img))==3)
        imshow(rgb2gray(img));
    else
        imshow(img);
    end
    alpha(0.2);
end
axis equal;
axis off;
hold all;
% set drawing options for the merging pairs
%cmap = brewermap(max(L),scheme);
if(L(i)<=n) % not yet clusterized
    c1 ='k';
else
    c1 = cmap(L(i),:);
end
if(L(j)<n) % not yet clusterized
    c2 ='k';
else
    c2 = cmap(L(j),:);
end
balli = L==L(i);
ballj = L==L(j);
% plot balls
if(~isempty(balli))
    scatter(X(1,balli),X(2,balli),mrkSizeClust,'o','MarkerFaceColor',c1,'LineWidth',eps,'MarkerFaceAlpha',.7,'MarkerEdgeAlpha',0);
end
if(~isempty(ballj))
    scatter(X(1,ballj),X(2,ballj),mrkSizeClust,'o','MarkerFaceColor',c2,'LineWidth',eps,'MarkerFaceAlpha',.7,'MarkerEdgeAlpha',0);
end
% plot not clusterized points
scatter(X(1,L<=n),X(2,L<=n),20,[0,0,0],'filled'); legend off;
% plot clusterized points
for u  = n+1:max(L)
    if(u ~=L(i) && u~=L(j))
       curSym = syms(mod(u-1,numel(syms))+1);
     scatter(X(1,L==u),X(2,L==u),mrkSize,cmap(u,:),curSym,'filled','MarkerEdgeColor',[0.2,0.2,0.2],'LineWidth',0.1,'MarkerFaceAlpha',0.8);
    end
end

% plot current clusters
scatter(X(1,L==L(i)),X(2,L==L(i)),mrkSizeBig,c1,'o','filled','MarkerEdgeColor','k','LineWidth',1.5,'MarkerFaceAlpha',1);
scatter(X(1,L==L(j)),X(2,L==L(j)),mrkSizeBig,c2,'o','filled','MarkerEdgeColor','k','LineWidth',1.5,'MarkerFaceAlpha',1);


end


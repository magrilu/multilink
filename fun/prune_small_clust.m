function [ C, card ] = prune_small_clust(C, theta )
%PRUNE_SMALL_CLUST prune small clusters: set to zero segments with
% theta points or less
 lbl = sort(unique(C));
 card = numel(lbl);
 
 for i = 1:numel(lbl)
    id = (C==lbl(i));
    card(i) = sum(id);
    if(card(i)<=theta)
        C(id)=0;
    end
        
 end
C(C~=0)=grp2idx(C(C~=0));
 
end


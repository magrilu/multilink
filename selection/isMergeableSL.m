function [ok] = isMergeableSL(P, L, i, j)
isInCi = L==L(i);
isInCj = L==L(j);
Pa = min(P(isInCi,:),[],1);
Pb = min(P(isInCj,:),[],1);

% compute jaccard index
if(all(Pa==0))
    jacc = 0;
elseif(all(Pb==0))
    jacc = 0;
else
    s = sum(Pa.*Pb);
    jacc = s/(Pa*Pa'+ Pb*Pb' - s);
end


ok = jacc> 0.0;

end
function [b] = is_invalid_flat(mss,delta)
%IS_INVALID_SUBSPACE 
b = rank(mss-mean(mss,2))<delta;

end


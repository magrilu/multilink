function F = fund_lin(m1,m2,w)
%FUND_LIN  Fundamental matrix with 8-points algorithm
    
    if(nargin <3)
        w = ones(size(m1,2),1);% weights
    end 
 
    % pre-conditioning
    [T1,m1] = precond(m1);  
    [T2,m2] = precond(m2);
     
    F = eight_points(m1, m2, w);
    
    % apply the inverse scaling
    F = T2' * F * T1;
    
    % enforce singularity of  F
    [U,D,V] = svd(F);
    D(3,3) = 0; F = U *D*V';
    
end














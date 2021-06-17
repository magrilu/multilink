function  y  = is_fundamental_degen( X, theta , cardmss )
%ISHDEGENERATE Check if a sample of 8 points X is H-degenerate w.r.t the
%fundamental matrices F. A sample is H_degenerate if there are 5 points in it
%related by an homography
%   INPUT:
%          X: 6x8 matrix-sample of eight points in homogeneous coordinate
%          in two views
%          y: boolean- true if X is H-degenerate.
%          theta: parameter vectors of fundamental matrices in columns
%   References: Two view geometry estimation unaffected by a dominant plane


%%

 tol = 1e-3;


if nargin==2
    cardmss=8;
end

N=size(theta,2);
y=false(1,N);
for I=1:N
    
    F = reshape(theta(:,I),3,3);
    %T = combntns(1:cardmss,5);
    T = nchoosek(1:cardmss,5);
    k = size(T,1);
   
    e = epipole(F');
    A  = star(e)*F;
    b=nan(3,1);
    
    
    
    
    for i = 1:k
        t = T(i,:);
        
        M = X(1:3, t(1:3))' ;
        if(rcond(M)<1e-10)
          y(I)=true;
          break;
        end
        for j = 1:3
            c = cross(X(4:6,t(j)),e);
            b(j) =  cross(X(4:6,t(j)),A*X(1:3,t(j)))'*c* ( norm( c ) )^(-2);
        end
        H = A - e*(M\b)'; % homography compatible with the triplet in t;
        % ceck for consistency with the other points
        res = nan(1,2); %true if point is consistent with the homography
        for j = 4:5
            %res(j-3) = norm(cross(H*X(1:3,t(j)),X(4:6,t(j))))
            a=H*X(1:3,t(j)); a=a./a(3);
            b=X(4:6,t(j));
            res(j-3) = norm(a-b);
        end
        
        y(I)=all(res<tol);
      
        if(y(I))
            break;
        end
    end
end
end



function a = fit_affinity(mss)
% return the parameter of the affinity matrix a = A(:);
% A is a 3x3 matrix 

[H A b] = affineLS(mss(1:2,:), mss(4:5,:));
a = H(:);

%%
% m1 = mss(1:2,:);
% m2 = mss(4:5,:);
% n = size(mss,2);
% AA = zeros(2*n,6); % coefficient matrix
% b = zeros(2*n,1);
% row = 1;
% for i = 1:n
%     AA(row,:)   = [m1(1,i),m1(2,i), 0,0, 1,0];
%     AA(row+1,:) = [0, 0, m1(1,i),m1(2,i), 0,1];
%     b(row) = m2(1,i);
%     b(row+1) = m2(2,i);
%     row = row+2;
% end
% 
% sol = AA\b;
% % arrange the solution in a 3x3 matrix
% A = [sol(1),sol(2),sol(5); sol(3),sol(4),sol(6); 0, 0 ,1];
% a = A(:);

end
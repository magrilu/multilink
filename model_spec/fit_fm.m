function f = fit_fm(mss)
%FIT_FM 
 %f = fund(mss(1:3,:),mss(4:6,:));
  f = fund_lin(mss(1:2,:),mss(4:5,:));
  f = f(:);
 %f = fm(mss(1:2,:),mss(4:5,:));
end


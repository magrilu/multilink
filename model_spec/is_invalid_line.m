function flg = is_invalid_line(h)
k = normc(h);
flg = (k(1)<1e-8 && k(2)<1e-8 ) || any(isnan(h));
end
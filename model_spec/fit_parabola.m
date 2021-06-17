function h = fit_parabola(X)
    h = polyfit(X(1,:),X(2,:),2);
    h = h(:);
end
function diffY_lags = predAR(diffY,phi,p,do_neg_SS)
    diffY_lags = diffY;
    for ip = 1:p
        diffY_lags((1+ip):end) = diffY_lags((1+ip):end) - phi(1+ip)*diffY(1:(end-ip));
    end
    if do_neg_SS
        diffY_lags = -sum(diffY_lags.^2);
    end
end
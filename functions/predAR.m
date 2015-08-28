function diffY_lags = predAR(diffY,phi,p,do_neg_SS)
    diffY_lags1 = diffY;
    for ip = 1:p
        diffY_lags1((1+ip):end) = diffY_lags1((1+ip):end) - phi(1+ip)*diffY(1:(end-ip));
    end
    if do_neg_SS
        diffY_lags1 = -sum(diffY_lags1.^2);
    end
% end

mask = triu(ones(p));

x = repmat(diffY(1:p),p,1).*mask;
% 
% diffY_lags = -sum((diffY-sum(bsxfun(@times,x,phi(2:p+1)'))).^2);

x = [sum(bsxfun(@times,x,phi(2:p+1)')) zeros(1,length(diffY)-p)];
sum(phi(2:p+1))
diffY_lags = -sum((diffY - sum(phi(2:p+1))*diffY + x).^2);



abs(diffY_lags - diffY_lags1)
function diffY_lags1 = predAR(diffY,phi,p,do_neg_SS,circ_ind,mask)
    diffY_lags1 = diffY;
    for ip = 1:p
        diffY_lags1((1+ip):end) = diffY_lags1((1+ip):end) - phi(1+ip)*diffY(1:(end-ip));
    end
    if do_neg_SS
        diffY_lags1 = -sum(diffY_lags1.^2);
    end
% end

% x = diffY(circ_ind).*mask;
% % 
% % diffY_lags = -sum((diffY-sum(bsxfun(@times,x,phi(2:p+1)'))).^2);
% 
% % Zx = [sum(bsxfun(@times,x,phi(2:p+1)')) zeros(1,length(diffY)-p)];
% % sum(phi(2:p+1))
% diffY_lags = -sum((diffY - sum(bsxfun(@times,x,phi(2:p+1)'))).^2);



% abs(diffY_lags - diffY_lags1)
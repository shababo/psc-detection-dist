function ef=genEfilt(tau,T)
% compute exponential filter(s)
%random stuff?
% D = 1;
% p1 = exp(D*-1/tau(1));
% p2 = exp(D*-1/tau(2));
% gamma_1 = p1+p2;
% gamma_2 = -p1*p2;

%% new way
ef_d = exp(-(0:T)/tau(2));

ef_h = -exp(-(0:T)/tau(1));

%compute maximum:
to = (tau(1)*tau(2))/(tau(2)-tau(1))*log(tau(2)/tau(1)); %time of maximum
max_val = exp(-to/tau(2))-exp(-to/tau(1)); %maximum
ef = {ef_h/max_val ef_d/max_val};


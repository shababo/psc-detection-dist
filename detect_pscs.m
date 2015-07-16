function detect_pscs(trace_file,param_file,param_ind,noise_type)

load(trace_file,'traces');
load(param_file,'a_min','p_spike','tau1_min','tau1_max','tau2_min','tau2_max');

param_dims = [length(a_min) length(p_spike) length(tau1_min) length(tau1_max) length(tau2_min) length(tau2_max)];
[a_min_i, p_spike_i, tau1_min_i, tau1_max_i, tau2_min_i, tau2_max_i] = ...
    ind2sub(param_dims,param_ind);

params.a_min = a_min(a_min_i);
params.p_spike = p_spike(p_spike_i);
params.tau1_min = tau1_min(tau1_min_i);
params.tau1_max = tau1_max(tau1_max_i);
params.tau2_min = tau2_min(tau2_min_i);
params.tau2_max = tau2_max(tau2_max_i);


params.dt = 1/20000;

% noise_types
gaussian = 1; line = 2; ar2 = 3;


% traces = traces_1_perm(2,start_t:end_t);

results = struct();

for trace_ind = 1:size(traces,1)
    
    trace_ind
    trace = max(traces(trace_ind,:)) - traces(trace_ind,:);

    % figure;plot(trace)

    % tGuess=[15 20];
    % tau = [3 9];
    %tGuess=[280 430 1345];
    %tGuess=[1345];
    tic
    tGuess = find_pscs_new(traces(trace_ind,:), dt, .002, 10, 0, 0);
    
    disp(['Starting events: ' num2str(length(tGuess))])
    
    tau = [5 35];
    switch noise_type
        case gaussian
            [results(trace_ind).trials, results(trace_ind).mcmc results(trace_ind).params]  = sampleParams(trace,tau,tGuess,dt);
        case line
            [results(trace_ind).trials, results(trace_ind).mcmc results(trace_ind).params]  = sampleParams_linenoise(trace,tau,tGuess,dt);
        case ar2
            [results(trace_ind).trials, results(trace_ind).mcmc results(trace_ind).params]  = sampleParams_ARnoise(trace,tau,tGuess,params);
    end
    results(trace_ind).runtime = toc;

% change tau min max and prior (and double check amplitudes and baseline
% limits
% amplitude threshold probably will help/is necessary.

end


%% minimum error sample

for trace_ind = 1:size(traces,1);

    trials = results(trace_ind).trials;
    times = results(trace_ind).times;
    mcmc = results(trace_ind).mcmc;
    trace = max(traces(trace_ind,:)) - traces(trace_ind,:);

    errP = zeros(1,length(trials.curves));
    for i = 1:length(trials.curves)
        errP(i) = sum((trials.curves{i}-trace).^2);
    end
    % figure;plot(errP)

    [results(trace_ind).min_err_ind results(trace_ind).min_err_ind] = min(errP);
    
end


save('0521_nice_l23_detectresults.mat','results')
%%
% plot MAP

trace_ind = 1;

trials = results(trace_ind).trials;
times = results(trace_ind).times;
mcmc = results(trace_ind).mcmc;
trace = max(traces(trace_ind,:)) - traces(trace_ind,:);

errP = zeros(1,length(trials.curves));
for i = 1:length(trials.curves)
    errP(i) = sum((trials.curves{i}-trace).^2);
end
% figure;plot(errP)

[me mi] = min(errP);



    
figure
plot(-trace,'c')
hold on
plot(-trials.curves{mi},'b','LineWidth',2)
hold on;
if fit_noise
    plot(-curve_no_line_noise+20,'b','LineWidth',2);
    curve_no_line_noise = remove_line_noise(trials.curves{mi},trials.phi{mi},trials.a_s{mi},60*dt,trace);
    hold on;
end
plot(-trace-20,'b')
hold off
axis tight
ylim([-50 50])
% axis off


%% plot samples
PlotMCMC

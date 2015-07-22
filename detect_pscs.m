function detect_pscs(trace_file,param_file,param_ind,noise_type)

rng(1234)

addpath(genpath('/vega/stats/users/bms2156/psc-detection'));

maxNumCompThreads(12)
matlabpool(12)

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

tic

parfor trace_ind = 1:size(traces,1)
    
    disp(['trace_ind = ' num2str(trace_ind)])
    trace = max(traces(trace_ind,:)) - traces(trace_ind,:);

    % figure;plot(trace)

    % tGuess=[15 20];
    % tau = [3 9];
    %tGuess=[280 430 1345];
    %tGuess=[1345];
%     tic
    tGuess = find_pscs(traces(trace_ind,:), params.dt, .002, 2, 1, 0, 0);
    
%     disp(['Starting events: ' num2str(length(tGuess))])
    
    tau = [5 35];
    switch noise_type
        case gaussian
            [results(trace_ind).trials, results(trace_ind).mcmc results(trace_ind).params]  = sampleParams(trace,tau,tGuess,params);
        case line
            [results(trace_ind).trials, results(trace_ind).mcmc results(trace_ind).params]  = sampleParams_linenoise(trace,tau,tGuess,params);
        case ar2
            [results(trace_ind).trials, results(trace_ind).mcmc results(trace_ind).params]  = sampleParams_ARnoise(trace,tau,tGuess,params);
    end
%     runtime = toc
%     results(trace_ind).runtime = runtime;    
    disp(['trace_ind = ' num2str(trace_ind) ', done in ' num2str(runtime) ' secs!'])

% change tau min max and prior (and double check amplitudes and baseline
% limits
% amplitude threshold probably will help/is necessary.

end

runtime = toc;
matlabpool close

%% minimum error sample

for trace_ind = 1:size(traces,1);

    [results(trace_ind).min_err, results(trace_ind).min_err_ind] = min(results(trace_ind).trials.obj);
    results(trace_ind).noise_type = noise_type;
    
end

savename = ['/vega/stats/users/bms2156/psc-detection/data/detection-results-' num2str(size(traces,1)) '-' num2str(noise_type) '-' num2str(param_ind) '.mat'];
save(savename,'results','runtime','noise_type','param_ind')


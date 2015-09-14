function detect_pscs(traces,params)

params = get_params;

if params.tau_min >= params.tau_max
    results = 'infeasible parameter set';
    save(params.savename,'results','params')
    return
end

if params.rand
    rng(params.seed)
end

addpath(genpath(params.source_path));

maxNumCompThreads(12)
matlabpool(12)

if ~isnumeric(traces)
    load(traces,'traces');
end

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


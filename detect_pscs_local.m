function detect_pscs_local(trace_file,params_in,param_ind,noise_type,event_sign)

rng(1234)

delete(gcp('nocreate'))
this_pool = parpool();



traces = [];
load(trace_file,'traces');


start_ind = 1;
duration = size(traces,2) - 1; %2150;
traces = event_sign*traces(:,start_ind:(start_ind + duration));

params.start_ind = start_ind; params.duration = duration;

% traces = traces(1,:);

if isnumeric(params_in)
    
    params.a_min = params_in(1);
    params.tau1_min = params_in(2);
    params.tau1_max = params_in(3);
    params.tau2_min = params_in(4);
    params.tau2_max = params_in(5);
    params.p_spike = params_in(6);
    
elseif isstr(params_in)
    
    load(params_in,'a_min','p_spike','tau_min','tau_max');


    param_dims = [length(a_min) length(p_spike) length(tau_min) length(tau_max)];
    [a_min_i, p_spike_i, tau_min_i, tau_max_i] = ...
        ind2sub(param_dims,param_ind);


    params.a_min = a_min(a_min_i);
    params.p_spike = p_spike(p_spike_i);
    params.tau1_min = tau1_min(tau1_min_i);
    params.tau1_max = tau1_max(tau1_max_i);
    params.tau2_min = tau2_min(tau2_min_i);
    params.tau2_max = tau2_max(tau2_max_i);
end


savename = ['data/local_test_' num2str(noise_type) '_' num2str(param_ind)  '.mat'];


if params.tau1_min >= params.tau1_max || params.tau2_min >= params.tau2_max
    results = 'infeasible parameter set';
    savename = [savename(1:end-4) '-z.mat'];
    save(savename,'results')
    return
end

params.dt = 1/20000;
params.p = 4;

% noise_types
gaussian = 1; line = 2; ar2 = 3;


% traces = traces_1_perm(2,start_t:end_t);

results = struct();
disp(size(traces,1));

% p = Par(size(traces,1));
% tic
parfor trace_ind = 1:size(traces,1)
%     
    disp(['trace_ind = ' num2str(trace_ind)])
    trace = max(traces(trace_ind,:)) - traces(trace_ind,:);

% tic
%     Par.tic
    tGuess = find_pscs(traces(trace_ind,:), params.dt, .002, params.a_min, 1);
    disp(['Starting events: ' num2str(length(tGuess))])
    
    tau = [mean([params.tau1_min params.tau1_max]) mean([params.tau2_min params.tau2_max])];
    switch noise_type
%         case gaussian
%             [results(trace_ind).trials, results(trace_ind).mcmc, results(trace_ind).params]  = sampleParams(trace,tau,tGuess,params);
%         case line
%             [results(trace_ind).trials, results(trace_ind).mcmc, results(trace_ind).params]  = sampleParams_linenoise(trace,tau,tGuess,params);
        case ar2
            [results(trace_ind).trials, results(trace_ind).mcmc, results(trace_ind).params]  = sampleParams_ARnoise_splittau(trace,tau,tGuess,params);
    end
% runtime = toc


%     p(trace_ind) = Par.toc;
%     results(trace_ind).runtime = p(trace_ind).ItStop - p(trace_ind).ItStart;

% change tau min max and prior (and double check amplitudes and baseline
% limits
% amplitude threshold probably will help/is necessary.

end
% stop(p)


delete(this_pool)

% for i = 1:length(results)
%     disp(results(i).runtime)
% end



%% map sample

for trace_ind = 1:size(traces,1);


    [results(trace_ind).min_err, results(trace_ind).min_err_ind] = min(results(trace_ind).trials.obj);
    
end

% savename = ['/vega/stats/users/bms2156/psc-detection/data/detection-results-' regexprep(mat2str(clock),'[| |\]|\d\d\.\d*','')];
save(savename,'results')





function infer_events(traces, params)

if params.rand == 1
    rng(params.seed)
end


delete(gcp('nocreate'))
this_pool = parpool();


traces = [];
load(trace_file,'traces');

if ~isfield(params,'start_ind')
    params.start_ind = 1;
end

if ~isfield(params,'duration')
    params.duration = size(traces,2);
end

traces = params.event_sign*traces(:,params.start_ind:(params.start_ind + params.duration - 1));


if isfield(params,'traces_ind')
    traces = traces(params.traces_ind,:);
end

if params.tau1_min >= params.tau1_max || params.tau2_min >= params.tau2_max
    results = 'infeasible parameter set';
    params.savename = [params.savename(1:end-4) '-z.mat'];
    save(params.savename,'results')
    return
end

results = struct();
disp(['About to run inference on :' num2str(size(traces,1)) ' traces...']);

parfor trace_ind = 1:size(traces,1)
%     
    disp(['Starting trace #' num2str(trace_ind)])
    trace = event_sign*traces(trace_ind,:);
    trace = trace - min(trace);

    event_times_init = template_matching(-1*params.event_sign*traces(trace_ind,:), params.dt,...
        params.init_method.tau, params.init_method, params.conv_thresh);
    
    tau = [mean([params.tau1_min params.tau1_max]) mean([params.tau2_min params.tau2_max])];

    [results(trace_ind).trials, results(trace_ind).mcmc, results(trace_ind).params]  = sampleParams_ARnoise_splittau(trace,tau,event_times_init,params);


end

delete(this_pool)

% map sample
for trace_ind = 1:size(traces,1);

    [results(trace_ind).min_err, results(trace_ind).min_err_ind] = min(results(trace_ind).trials.obj);
    
end

save(params.full_save_string,'results')





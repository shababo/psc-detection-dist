function infer_events(params)

if params.rand == 1
    rng(params.seed)
end

if params.cluster
    addpath('/vega/stats/users/bms2156/psc-detection/functions')
end


if ~isfield(params,'start_ind')
    params.start_ind = 1;
end

try
    load(params.traces_filename,'traces')
catch e
    params.traces_filename = 'data/for-paper/direct-stim-w-events-real.mat';
    load('data/for-paper/direct-stim-w-events-real.mat')
end



% the file were your traces are, traces should be saved in this mat file in
% an N x T matrix called 'traces' where N = number of traces and T = number
% of samples
% load(params.traces_filename,'traces');

if params.is_grid
    [traces, rebuild_map] = stack_traces(traces);
    params.rebuild_map = rebuild_map;
end

% assignin('base','rebuild_map',rebuild_map)
% return

if ~isfield(params,'duration')
    params.duration = size(traces,2);
end

% grab section in time
traces = traces(:,params.start_ind:(params.start_ind + params.duration - 1));

% grab subset of traces
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
disp(['About to run inference on: ' num2str(size(traces,1)) ' traces...']);

if params.par
    
    if params.cluster

        addpath(genpath(params.source_path));

        maxNumCompThreads(12)
        matlabpool(12)

    else

        delete(gcp('nocreate'))
        this_pool = parpool();

    end
    
    load_struct = load(params.init_method.template_file);
    template = load_struct.template;
    
    parfor trace_ind = 1:size(traces,1)
    %     
        disp(['Starting trace #' num2str(trace_ind)])
        trace = params.event_sign*traces(trace_ind,:);
        trace = trace - min(trace);

%         event_times_init_old = template_matching(-1*params.event_sign*traces(trace_ind,:), params.dt,...
%             params.init_method.tau, params.init_method.amp_thresh, params.init_method.conv_thresh);
        
        
        
        nfft = length(trace) + length(template);
        [filtered_trace, event_times_init,event_sizes_init] = wiener_filter(trace,template,params.init_method.ar_noise_params,...
            nfft, params.dt, params.init_method.theshold, params.init_method.min_interval);
        event_times_init
        event_sizes_init
%         assignin('base','event_times_init_old',event_times_init_old)
%         assignin('base','event_times_init',event_times_init)
        results(trace_ind).event_times_init = event_times_init;
        results(trace_ind).filtered_trace = filtered_trace;
        results(trace_ind).event_sizes_init = event_sizes_init;
        
        tau = [mean([params.tau1_min params.tau1_max]) mean([params.tau2_min params.tau2_max])]/params.dt;
        
        if isfield(params,'init_only') && ~params.init_only
            if params.direct_stim
    %             event_times_init = ceil(length(trace)*rand(1,length(trace)*params.p_spike));
                [results(trace_ind).trials, results(trace_ind).mcmc]  = sampleParams_ar_2taus_directstim(trace,tau,event_times_init,params);
            else
    %             event_times_init = template_matching(-1*params.event_sign*traces(trace_ind,:), params.dt,...
    %                 params.init_method.tau, params.init_method.amp_thresh, params.init_method.conv_thresh);
                [results(trace_ind).trials, results(trace_ind).mcmc]  = sampleParams_ARnoise_splittau(trace,tau,event_times_init,params);
            end
        end
    end
    
    
    if params.cluster
        matlabpool close
    else
        delete(this_pool)
    end

else
    
            
    load_struct = load(params.init_method.template_file);
    template = load_struct.template;

    for trace_ind = 1:size(traces,1)
        disp(['Starting trace #' num2str(trace_ind)])
        trace = params.event_sign*traces(trace_ind,:);
        trace = trace - min(trace);

%         event_times_init_old = template_matching(-1*params.event_sign*traces(trace_ind,:), params.dt,...
%             params.init_method.tau, params.init_method.amp_thresh, params.init_method.conv_thresh);


        
        nfft = length(trace) + length(template);
        [~, event_times_init,event_sizes_init] = wiener_filter(params.event_sign*trace,template,params.init_method.ar_noise_params,...
            nfft, params.dt, params.init_method.theshold, params.init_method.min_interval);
        event_times_init
        event_sizes_init
        results(trace_ind).event_times_init = event_times_init;
%         results(trace_ind).filtered_trace = filtered_trace;
        results(trace_ind).event_sizes_init = event_sizes_init;
%         assignin('base','event_times_init_old',event_times_init_old)
%         assignin('base','event_times_init',event_times_init)
        
        tau = [mean([params.tau1_min params.tau1_max]) mean([params.tau2_min params.tau2_max])]/params.dt;
        
        if isfield(params,'init_only') && ~params.init_only

            if params.direct_stim
    %             event_times_init = ceil(length(trace)*rand(1,length(trace)*params.p_spike));
                [results(trace_ind).trials, results(trace_ind).mcmc]  = sampleParams_ar_2taus_directstim(trace,tau,event_times_init,params);
            else
                [results(trace_ind).trials, results(trace_ind).mcmc]  = sampleParams_ARnoise_splittau(trace,tau,event_times_init,params);
            end
        end
    end
end

if isfield(params,'init_only') && ~params.init_only
    disp('finding min err...')
    % map sample
    for trace_ind = 1:size(traces,1)

        [results(trace_ind).map, results(trace_ind).map_ind] = max(results(trace_ind).trials.obj);

    end
end

% if params.is_grid
%     results_grid = unstack_results(results, rebuild_map);
% end

% results = results_grid;

disp('saving...')
save(params.full_save_string,'results','params')

disp('done')



function run_posterior_sampler(params)

% if we want to seed the rng before running
if params.rand == 1
    rng(params.seed)
end

% add source if it's not already on path - useful for running on a cluster
if params.do_addpath
    addpath(genpath(params.source_path));
end

% attempt to load data
try
    load(params.traces_filename,'traces')
catch e
    disp('Could not load traces file, aborting...');
    return
end

% initialize results struct
results = struct();

% output progress to user
disp(['About to run inference on: ' num2str(size(traces,1)) ' traces...']);

% load template for init method
load_struct = load(params.init_method.template_file);
template = load_struct.template;

% if we want to use parfor
if params.par
    
    % determine version year for MATLAB
    version_year = version('-release');
    version_year = str2double(version_year(1:4));
    
    % init parallel resources based on version
    if version_year <= 2013
        maxNumCompThreads(12)
        matlabpool(12)
    else
        delete(gcp('nocreate'))
        this_pool = parpool();
    end
    
    
    % run parfor over traces
    parfor trace_ind = 1:size(traces,1)
        
        % display progress
        disp(['Starting trace #' num2str(trace_ind)])
        
        % get trace, make events go up, make trace nonnegative
        trace = params.event_sign*traces(trace_ind,:);
        trace = trace - min(trace);

        % run wiener filter to init event times
        nfft = length(trace) + length(template);
        
        % run wiener filter
        [filtered_trace, event_times_init, event_sizes_init] = ...
            wiener_filter(trace,template,params.init_method.ar_noise_params,...
            nfft, params.dt, params.init_method.theshold, params.init_method.min_interval);
        
        % store initialization results
        results(trace_ind).event_times_init = event_times_init;
        results(trace_ind).filtered_trace = filtered_trace;
        results(trace_ind).event_sizes_init = event_sizes_init;
        
        % make init for taus
        tau = [mean([params.tau1_min params.tau1_max]) ...
               mean([params.tau2_min params.tau2_max])]/params.dt;
        taus_init = repmat(tau,length(event_times_init),1);
        
        % run sampler
        if ~params.init_only
            [results(trace_ind).trials, results(trace_ind).mcmc]  = ...
                sample_params(trace, params, event_times_init, ...
                event_sizes_init, taus_init);
        end
        
    end
    
    % shut down parallel resources
    if version_year <= 2013
        matlabpool close
    else
        delete(this_pool)
    end

% sample params for each trace serially
else

    for trace_ind = 1:size(traces,1)
        
        % display progress
        disp(['Starting trace #' num2str(trace_ind)])
        
        % get trace, make events go up, make trace nonnegative
        trace = params.event_sign*traces(trace_ind,:);
        trace = trace - min(trace);

        % run wiener filter to init event times
        nfft = length(trace) + length(template);
        
        [filtered_trace, event_times_init, event_sizes_init] = ...
            wiener_filter(params.event_sign*trace,template,params.init_method.ar_noise_params,...
            nfft, params.dt, params.init_method.theshold, params.init_method.min_interval);

        % store initialization results
        results(trace_ind).event_times_init = event_times_init;
        results(trace_ind).filtered_trace = filtered_trace;
        results(trace_ind).event_sizes_init = event_sizes_init;
        
        % create tau init
        tau = [mean([params.tau1_min params.tau1_max]) ...
               mean([params.tau2_min params.tau2_max])]/params.dt;
        taus_init = repmat(tau,length(event_times_init),1);
        
        % run sampler
        if ~params.init_only
            [results(trace_ind).trials, results(trace_ind).mcmc]  = ...
                sample_params(trace, params, event_times_init, ...
                event_sizes_init, taus_init);
        end
    end
end

% for each trace, find the map sample
if isfield(params,'init_only') && ~params.init_only
    disp('finding map sample...')
    
    % map sample
    for trace_ind = 1:size(traces,1)
        [results(trace_ind).map, results(trace_ind).map_ind] = max(results(trace_ind).trials.obj);
    end
    
end

% save it
disp('saving...')
save(params.full_save_string,'results','params')

% ...and we're out.
disp('done')



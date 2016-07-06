function params = get_params(base_params_path)


% load a params struct from a file to start with
if ~isempty(base_params_path)
    load(base_params_path);
    
% or create a new struct
else
    params = struct();
end

% are we using a cluster?
if ~isfield(params,'cluster')
    params.cluster = 1;
end

% do we need to put the source into the path for MATLAB?
if ~isfield(params,'do_addpath')
    params.do_addpath = params.cluster; % generally we use this when we also are using the cluster
end

% if so, what is the source path?
if ~isfield(params,'source_path')
    params.source_path = '/vega/stats/users/bms2156/psc-detection';
end

% use MATLAB's parfor across traces?
if ~isfield(params,'par')
    params.par = 1;
end

% only run initialization algorithm (usually Wiener Filter)
if ~isfield(params,'init_only')
    params.init_only = 0;
end

%% use an rng seed

% seed the rng?
if ~isfield(params,'rand')
    params.rand = 1;
end

% pick a seed
if ~isfield(params,'seed')
    params.seed = 12341;
end


%% data params

% time in seconds per sample
if ~isfield(params,'dt')
    params.dt = 1/20000;
end

% direction/sign of events: upward is 1 (e.g. ipscs, ca imaging), downard is -1
% (e.g. epscs)
if ~isfield(params,'event_sign')
    params.event_sign = 1;
end


%% inference params

% event amplitude bounds - a_min needs to be less than a_max, but a_min
% does not have to be greater than or equal to 0. That is,
% you can simultaneously detect upward and downward events.
if ~isfield(params,'a_max')
    params.a_max = Inf;
end
if ~isfield(params,'a_min')
    params.a_min = 5;
end

% baseline bounds - i.e. the holding current
if ~isfield(params,'b_min')
    params.b_min = -200;
end
if ~isfield(params,'b_max')
    params.b_max = 200;
end

% min and max for "rise time" in seconds
if ~isfield(params,'tau1_min')
    params.tau1_min = 10/20000;
end
if ~isfield(params,'tau1_max')
    params.tau1_max = 100/20000;
end

% min and max for "decay time" in seconds
if ~isfield(params,'tau2_min')
    params.tau2_min = 100/20000;
end
if ~isfield(params,'tau2_max')
    params.tau2_max = 700/20000;
end

% how long to make kernel in samples
if ~isfield(params,'event_samples')
    params.event_samples = 6*params.tau2_max/params.dt;
end

% poisson rate in events/sample
if ~isfield(params,'p_event')
    params.p_event = 1e-4;
end


% ar noise model parameters
% p, the number of timesteps for filter
if ~isfield(params,'p')
    params.p = 2; % how many time steps to regress on
end
% mean for MVN prior on phi - we find that the ar parameters are - to some extent - 
% headstage dependent for voltage clamp recordings
if ~isfield(params,'phi_0')
    params.phi_0 = [1.0, -0.30]';
end
% inverse covariance/precision matrix for MVN prior on phi
if ~isfield(params,'Phi_0')
    params.Phi_0 = 10*eye(params.p); %inverse covariance
end
% prior parameters for inverse gamma on sigma_sq
if ~isfield(params,'nu_0')
    params.nu_0 = 0;
end
if ~isfield(params,'sig2_0')
    params.sig2_0 = .1;
end
% initial value for variance for white noise input to filter
if ~isfield(params,'noise_var_init')
    params.noise_var_init = 3.5;
end

% it can be useful to estimate the noise from a small amount of data and
% then take it as known for further analysis within that recordings - or
% possibly even more generally across recordings. if the noise model is
% known, then enter it here.
if ~isfield(params, 'noise_known')
    params.noise_known = 1;
    if params.noise_known
        params.phi_known = [1.000000000000000 1.0 -.30];
        params.noise_var_known = 3.5;
    end
end

% select an subset of each trace, in samples, e.g., 1:1000,
% to use for noise estimateion. this can be useful when there are sections
% of the trace with a high-rate of events due to some stimulus
if ~isfield(params,'noise_est_subset')
    params.noise_est_subset = [];
end


%% sampling params


% how long to run the sampler
if ~isfield(params,'num_sweeps')
    params.num_sweeps = 2000;
end

% nubmer of burn in sweeps to run
if ~isfield(params,'burn_in_sweeps')
    params.burn_in_sweeps = 0;
end

% sampling event times proposal variance in seconds
if ~isfield(params,'time_proposal_var')
    params.time_proposal_var = 7.5e-04;
end

% rise time proposal variance in seconds
if ~isfield(params,'tau1_prop_std')
    params.tau1_prop_std = 2/20000;
end

% fall time proposal variance in seconds
if ~isfield(params,'tau2_prop_std')
    params.tau2_prop_std = 20/20000;
end

% amplitude proposal variance in pA or trace units
if ~isfield(params,'amp_prop_std')
    params.amp_prop_std = 3;
end

% baseline proposal variance in pA or trace units
if ~isfield(params,'baseline_prop_std')
    params.baseline_prop_std = 2;
end

% for each parameter we run an inner loop and collect extra samples - this
% can help balance the amount of samples we run for different parameters.
% only the last sample for each parameter is kept such that there this is
% still only one sample per sweep

% number of add/drop subsweeps per sweep
if ~isfield(params,'add_drop_sweeps')
    params.add_drop_sweeps = 10;
end

% number of event time subsweeps per sweep
if ~isfield(params,'time_sweeps')
    params.event_time_sweeps = 10;
end

% number of amplitude subsweeps per sweep
if ~isfield(params,'amp_sweeps')
    params.amp_sweeps = 5;
end

% number of baseline/holding current subsweeps per sweep
if ~isfield(params,'baseline_sweeps')
    params.baseline_sweeps = 1;
end

% number of rise time subsweeps per sweep
if ~isfield(params,'tau1_sweeps')
    params.tau1_sweeps = 1;
end
% number of fall time subsweeps per sweep
if ~isfield(params,'tau2_sweeps')
    params.tau2_sweeps = 1;
end

% number of samples for window around new event time to guide amplitude
% initialization - essentially the max or min within that window
if ~isfield(params,'a_init_window')
    params.a_init_window = 50;
end

% do not detect events within this many seconds of another event
if ~isfield(params,'exclusion_bound')
    params.exclusion_bound = 10/20000;
end

%% initialization method - best choice is usually a wiener filter

if ~isfield(params,'init_method')

    % ipsc template file
    params.init_method.template_file = 'data/ipsc-template.mat';
    % ipsc
%     params.init_method.template_file = 'data/epsc-template.mat';

    % noise parameters for Wiener Filter
    params.init_method.ar_noise_params.sigma_sq = 3.0;
    params.init_method.ar_noise_params.phi = [1.000000000000000, 1.0, -0.20];
    
    % threshold output of filter
    params.init_method.theshold = 2.25; % in std devs of filtered trace
    params.init_method.min_interval = 20; % in samples
end


%% paths for data
if ~isfield(params,'traces_filename')
    params.traces_filename = ...
        ['data/sample_data.mat'];
end

if ~isfield(params,'savepath')
    params.savepath = '';
end

% use the four digit identifier to run different versions of inference on
% the same data. if the name is already taken, the algorithm will abort
params.savename = [params.traces_filename(1:end-4) '-0000.mat'];
params.full_save_string = [params.savename];



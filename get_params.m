function params = get_params
rng(1234)

%% data params
% time in seconds per sample
params.dt = 1/20000;

% direction/sign of events: upward is 1 (e.g. ipscs, ca imaging), downard is -1
% (e.g. epscs)
params.event_sign = -1;

%% subtraces

% first sample, if you want to start at 1, omit
% params.start_ind = .33*20000;

% if you want to go to the end of the traces, omit
% params.duration = 20000*.1;

% if you want all traces, omit
params.traces_ind = randsample(80,18);

%% inference params

% event amplitude bounds
params.a_max = Inf;
params.a_min = 2.5;

% baseline bounds
params.b_min = -50;
params.b_max = 50;

% event kernel params
params.feature_names = {'amplitude','tau 1','tau 2','time'};
% params.kernel = @kernel_function; ignore this
% min and max for "rise time" in seconds
params.tau1_min = 1/20000;
% params.tau1_max = 60/20000;
% params.tau2_min = 75/20000;
% params.tau2_max = 300/20000;
% params.event_samples = 6*params.tau2_max/params.dt;
% 
% % poisson/rate
% params.p_spike = 1e-3;
params.tau1_max = 30/20000;
% min and max for "decay time" in seconds
params.tau2_min = 10/20000;
params.tau2_max = 300/20000;
% how long to make kernel in samples
params.event_samples = 6*params.tau2_max/params.dt;

% poisson/rate - that is the probability of seeing a spike/sample
params.p_spike = 1e-9;



% ar noise model
params.p = 4; % how many time steps to regress on
params.phi_0 = zeros(params.p,1);
params.Phi_0 = 10*eye(params.p); %inverse covariance 3

params.noise_var_init = 5;

%% sampling params


% how long to run the sampler
params.num_sweeps = 300;
params.burn_in_sweeps = 0;

% sampling spike times
params.time_proposal_var = 10;

params.tau1_prop_std = 2/20000;
params.tau2_prop_std = 20/20000;

params.amp_prop_std = .2;

params.baseline_prop_std = 2;

params.add_drop_sweeps = 5;
params.spike_time_sweeps = 3;
params.amp_sweeps = 5;
params.baseline_sweeps = 1;
params.tau1_sweeps = 1;
params.tau2_sweeps = 1;

params.exclusion_bound = 1;
params.Dt = 1;
params.A = 1;
% params.b
%% template-matching initialization method
params.init_method.tau = .002; % min seconds
params.init_method.amp_thresh = 5;
params.init_method.conv_thresh = 1;


%% sourcefile (for cluster)

params.source_path = '/vega/stats/users/bms2156/psc-detection';


%% use an rng seed

params.rand = 1;
params.seed = 1234;



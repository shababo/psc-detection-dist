function params = get_params

%% savename

params.source_path = '/vega/stats/users/bms2156/psc-detection';
params.savepath = 
params.savename = 

%% use an rng seed

params.rand = 1;
params.seed = 1234;

%% subtraces

% currently using samples, but could use time...
params.start_ind = 1;

% if you want to go to the end of the traces, omit
% params.duration = 1000;

% if you want all traces, omit
% params.traces_ind = 1:10;

%% inference params

% amplitude
params.inference.a_max = Inf;
params.inference.a_min = 1;

% event kernel
params.kernel
params.inference.tau1_min = 1;
params.inference.tau1_max = 20;
params.inference.tau2_min = 20;
params.inference.tau2_max = 75;

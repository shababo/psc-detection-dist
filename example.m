% This code generates some fake data and runs the sampler on it

% for testing
rng(12341)

% number of traces to generate
K = 1;

% set parameters for generating data
params.T = 20000;
params.dt = 1/20000;
params.phi = [.80, -.12];
params.sigma_sq = 3.0;
params.baseline_bounds = [-100 0];
params.tau_r_bounds = [0.25 1.0]*1e-3;
params.tau_f_bounds = [0.0025    0.0075];
params.a_bounds = [.5 20];
params.rate = 20;
params.event_sign = -1;
params.event_direction = params.event_sign;

[traces, true_signal] = sample_traces(K,params);

% plot the noisy data
figure
plot_trace_stack(traces,40,zeros(size(traces,1),3),'-')
title('Noisy Data')

% plot the noise-free data
figure
plot_trace_stack(true_signal.traces,40,zeros(size(traces,1),3),'-')
title('Noise-free Data')

% save data to runresults
filename = 'data/test-data.mat';
save(filename,'traces','true_signal')

% build a parameter struct that will point to this file
params.traces_filename = filename;
% fill rest of struct with defaults
params = get_params(params);

% check that you're not writing over a previous results file
if exist(params.full_save_string, 'file') == 2
    disp('****ABORTING: THE REQUESTED RESULTS FILE NAME ALREADY EXISTS****')
    disp(params.full_save_string)
    return
end

% run sampler
run_posterior_sampler(params);




function params = get_params(varargin)

if ~isempty(varargin)
    params = varargin{1};
end

if ~isfield(params,'cluster')

    params.cluster = 0;

end

if ~isfield(params,'par')
    params.par = 1;
end
%% use an rng seed

if ~isfield(params,'rand')
    params.rand = 1;
end

if ~isfield(params,'seed')
    params.seed = 12341;
end

% rng(params.seed)

%%


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

%% subtraces

% first sample, if you want to start at 1, omit
if ~isfield(params,'start_ind')
%     params.start_ind = .3*20000;
end
% if you want to go to the end of the traces, omit
if ~isfield(params,'duration')
%     params.duration = 20000*.15;
end

% if you want all traces, omit
if ~isfield(params,'traces_ind')
%     params.traces_ind = randsample(80,18);

%     params.traces_ind = 1278;
%     params.traces_ind = 1;

%     params.traces_ind = 1:6;
end
%% inference params

% event amplitude bounds
if ~isfield(params,'a_max')
    params.a_max = Inf;
end
if ~isfield(params,'a_min')
    params.a_min = 2;

end

% baseline bounds
if ~isfield(params,'b_min')
    params.b_min = -50;
end
if ~isfield(params,'b_max')
    params.b_max = 50;
end

% event kernel params
if ~isfield(params,'feature_names')
    params.feature_names = {'amplitude','tau 1','tau 2','time'};
end
% params.kernel = @kernel_function; ignore this
% min and max for "rise time" in seconds
if ~isfield(params,'tau1_min')
    params.tau1_min = 5/20000;
end
% params.tau1_max = 60/20000;
% params.tau2_min = 75/20000;
% params.tau2_max = 300/20000;
% params.event_samples = 6*params.tau2_max/params.dt;
% 
% % poisson/rate
% params.p_spike = 1e-3;
if ~isfield(params,'tau1_max')
    params.tau1_max = 20/20000;
end
% min and max for "decay time" in seconds
if ~isfield(params,'tau2_min')
    params.tau2_min = 20/20000;
end
if ~isfield(params,'tau2_max')
    params.tau2_max = 150/20000;
end
% how long to make kernel in samples
if ~isfield(params,'event_samples')
    params.event_samples = 6*params.tau2_max/params.dt;
end

% poisson/rate - that is the probability of seeing a spike/sample
if ~isfield(params,'p_spike')
    params.p_spike = 1e-4;%1e-4;

end



% ar noise model
if ~isfield(params,'p')
    params.p = 2; % how many time steps to regress on
end
if ~isfield(params,'phi_0')
    params.phi_0 = zeros(params.p,1);
end
if ~isfield(params,'Phi_0')
    params.Phi_0 = 10*eye(params.p); %inverse covariance 3
end

if ~isfield(params,'noise_var_init')
    params.noise_var_init = 1.0;
end

if ~isfield(params,'noise_est_subset')
%     params.noise_est_subset = 1:5000;
%     params.noise_est_subset = 1:1000;
end

%% direct stim

if ~isfield(params,'direct_stim')
    params.direct_stim = 0;
end

if ~isfield(params,'stim_tau_rise')
    params.stim_tau_rise = .003;
end
if ~isfield(params,'stim_tau_fall')
    params.stim_tau_fall = .050;
end

if ~isfield(params,'stim_amp_std')
    params.stim_amp_std = 10; %pA
end

if ~isfield(params,'stim_amp_min')
    params.stim_amp_min = 0;
end

if ~isfield(params,'stim_amp_max')
    params.stim_amp_max = Inf;
end

if ~isfield(params,'stim_in')
    params.stim_in = [zeros(1,20*5) ones(1,20*20) zeros(1,1000-25*5)];
end

if ~isfield(params,'stim_amp_init')
    params.stim_amp_init = NaN;
end

if ~isfield(params,'stim_tau_rise_min')
    params.stim_tau_rise_min = .0005;
end

if ~isfield(params,'stim_tau_rise_max')
    params.stim_tau_rise_max = .010;
end

if ~isfield(params,'stim_tau_fall_min')
    params.stim_tau_fall_min = .005;
end

if ~isfield(params,'stim_tau_fall_max')
    params.stim_tau_fall_max = .200;
end

if ~isfield(params,'stim_tau_rise_std')
    params.stim_tau_rise_std = .0005;
end

if ~isfield(params,'stim_tau_fall_std')
    params.stim_tau_fall_std = .005;
end

if ~isfield(params,'stim_shape')
%     load('data/for-paper/chr2-stim-response.mat');
%     params.stim_shape = chr2_response;
    
    load('data/2P-Chrimson-10ms-template.mat');
    params.stim_shape = -template_clean;
    
%     params.stim_shape = [];
end


%% sampling params


% how long to run the sampler
if ~isfield(params,'num_sweeps')
    params.num_sweeps = 2000;
end
if ~isfield(params,'burn_in_sweeps')
    params.burn_in_sweeps = 0;
end

% sampling spike times
if ~isfield(params,'time_proposal_var')
    params.time_proposal_var = 15;
end

if ~isfield(params,'tau1_prop_std')
    params.tau1_prop_std = 2/20000;
end
if ~isfield(params,'tau2_prop_std')
    params.tau2_prop_std = 20/20000;
end

if ~isfield(params,'amp_prop_std')
    params.amp_prop_std = .3;
end
if ~isfield(params,'baseline_prop_std')
    params.baseline_prop_std = 2;
end

if ~isfield(params,'add_drop_sweeps')
    params.add_drop_sweeps = 10;
end
if ~isfield(params,'time_sweeps')
    params.spike_time_sweeps = 10;
end
if ~isfield(params,'amp_sweeps')
    params.amp_sweeps = 5;
end
if ~isfield(params,'baseline_sweeps')
    params.baseline_sweeps = 1;
end
if ~isfield(params,'tau1_sweeps')
    params.tau1_sweeps = 1;
end
if ~isfield(params,'tau2_sweeps')
    params.tau2_sweeps = 1;
end

if ~isfield(params,'exclusion_bound')
    params.exclusion_bound = 10;
end
if ~isfield(params,'Dt')
    params.Dt = 1;
end
if ~isfield(params,'A')
    params.A = 1;
end
% params.b
%% template-matching initialization method
if ~isfield(params,'init_method')
    params.init_method.tau = .002; % min seconds
    params.init_method.amp_thresh = 5;
    params.init_method.conv_thresh = 1;
end


%% sourcefile (for cluster)
if ~isfield(params,'source_path')
    params.source_path = '/vega/stats/users/bms2156/psc-detection';
end

%% filenames
if ~isfield(params,'traces_filename')
    if params.cluster

        params.traces_filename = '/vega/stats/users/bms2156/psc-detection/data/for-paper/all-evoked-ipscs.mat';

    else
        params.traces_filename = ...
            ['data/for-paper/doublet-data-01.mat'];

    end
end

if ~isfield(params,'savepath')
    if params.cluster
        params.savepath = '/vega/stats/users/bms2156/psc-detection/data';
    else
        params.savepath = 'data/';
    end
end
if ~isfield(params,'savename')
    if params.cluster
%         savefile_basename = '/simulated-epscs-1027-results-0000-pspike-%0.0e-amin-%0.0e-num_sweeps-%0.0e.mat';
%         params.savename = sprintf(savefile_basename,params.p_spike,params.a_min,params.num_sweeps);
%         params.savename = strrep(params.savename,'+','');
        params.savename = 'all-evoked-ipscs-0000.mat';
    else
        params.savename = 'doublet-data-01-0000.mat';
    end
end

params.full_save_string = [params.savepath '/' params.savename];



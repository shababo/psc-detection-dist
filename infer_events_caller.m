function infer_events_caller(varargin)
% if you are loading params from a mat file (for example a previous run)
% if you are going to use the get_params function comment out
% paramfile = 'data/emx_ipscs_3cells3traces-results-0008.mat';

% varargin{1}: paramfile

params = struct();
        
% load params from file if given
if ~isempty(varargin) && ~isempty(varargin{1})
    varargin{1}
    load(varargin{1},'params');
%     [pathname, filename] = fileparts(params.traces_filename);
%     pathname = '/media/shababo/Layover/projects/mapping/code/psc-detection/data/for-paper/';
%     params.traces_filename = [pathname filename];

    params.savename = 'evoked_pscs-new-results-0011.mat';
% %     params.stim_in = [zeros(1,20*10) ones(1,20*10) zeros(1,2000-20*20)];%linspace(0,1,20*10)
%     % ar noise model
%     params.p = 3; % how many time steps to regress on
%     params.phi_0 = zeros(params.p,1);
%     params.Phi_0 = 10*eye(params.p); %inverse covariance 3
    params.a_min = .5;
%     load('data/for-paper/chr2-stim-response.mat');
%     params.stim_shape = chr2_response;
%     params.par = 0;

    params.p_spike = 1e-3;
%     params.noise_est_subset = 1:800;
%     params.a_min = 1;
%     params.start_ind = 7200;
%     params.duration = 800;
%     params.tau1_max = 25;
%     params.tau2_min = 5;
%     params.tau2_max = 200;
    
end

% fill with default params
params = get_params(params);

% check that you're not writing over a previous results file
if exist(params.full_save_string, 'file') == 2
    disp('****ABORTING: THE REQUESTED RESULTS FILE NAME ALREADY EXISTS****')
    disp(params.full_save_string)
    return
end

% test savefile before we do inference
results = 'temp'; 
save(params.full_save_string,'results','params','-v7.3')

% infer the events!!!
infer_events(params)

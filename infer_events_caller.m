function infer_events_caller

% 
% load('data/emx_ipscs_3cells3traces-results-0008.mat','params');
% params.num_sweeps = 500;

params = get_params;

params.traces_filename = 'data/pv_ipscs_3cells3traces.mat';

load(params.traces_filename,'traces');

params.savepath = 'data/';
params.savename = 'pv_ipscs_3cells3traces-results-0012.mat';
params.full_save_string = [params.savepath '/' params.savename];

if exist(params.full_save_string, 'file') == 2
    disp('****ABORTING: THE REQUESTED RESULTS FILE NAME ALREADY EXISTS****')
    return
end

infer_events(traces,params)
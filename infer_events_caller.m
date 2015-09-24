function infer_events_caller

% if you are loading params from a mat file (for example a previous run)
% if you are going to use the get_params function comment out
paramfile = 'data/emx_ipscs_3cells3traces-results-0008.mat';

if exist(paramfile,'var')
    load(paramfile,'params');
else
    params = get_params;
end

% the file were your traces are, traces should be saved in this mat file in
% an N x T matrix called 'traces' where N = number of traces and T = number
% of samples
params.traces_filename = 'data/vip_ipscs_3cells3traces.mat';
load(params.traces_filename,'traces');

% the directory to save the results in
params.savepath = 'data/';
% name of the file to save results in
params.savename = 'vip_ipscs_3cells3traces-results-0000.mat';
params.full_save_string = [params.savepath '/' params.savename];

% check that you're not writing over a previous results file
if exist(params.full_save_string, 'file') == 2
    disp('****ABORTING: THE REQUESTED RESULTS FILE NAME ALREADY EXISTS****')
    return
end

% infer the events!!!
infer_events(traces,params)
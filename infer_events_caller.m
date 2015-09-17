function infer_events_caller


load('data/evoked-pscs-strong-results-0000.mat','params');
params.num_sweeps = 4000;
% params = get_params;
% params.traces_filename = 'data/evoked_pscs.mat';

load(params.traces_filename,'traces');

params.savepath = 'data/';
params.savename = 'evoked-pscs-strong-results-0001.mat';
params.full_save_string = [params.savepath '/' params.savename];

params.num_sweeps = 4000;

if exist(params.full_save_string, 'file') == 2
    disp(['****ABORTING: THE REQUESTED RESULTS FILE NAME ALREADY EXISTS****'])
    return
end

infer_events(traces,params)
function infer_events_caller(varargin)
% if you are loading params from a mat file (for example a previous run)
% if you are going to use the get_params function comment out
% paramfile = 'data/emx_ipscs_3cells3traces-results-0008.mat';

% varargin{1}: paramfile

params = struct();

% load params from file if given
if ~isempty(varargin) && ~isempty(varargin{1})
    load(varargin{1},'params');
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
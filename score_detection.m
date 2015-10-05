function score_detection(resultsfile, tolerance, plot_results)

load(resultsfile)
params.timing_score_tolerance = tolerance;
% load(params.traces_filename,'true_event_times')
load('/home/shababo/projects/mapping/code/psc-detection/data/simulated-data-longer-traces-epsc.mat')

num_traces = length(results);

for i = 1:num_traces
    
    map_ind = results(i).map_ind;
    est_times = results(i).trials.times{map_ind};
    [timing_score(i).correct_inds, timing_score(i).correct_err] = ...
        score_detection_trace(est_times, true_event_times{i}, tolerance);
    
end

save(resultsfile,'timing_score','-append')
    
if plot_results
    
    figure
    plot_score_detection(results,true_event_times,timing_score)
    
end
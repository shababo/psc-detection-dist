function score_detection_deconv(resultsfile, tolerance)

load(resultsfile)

% load(params.traces_filename,'true_event_times')
load('/home/shababo/projects/mapping/code/psc-detection/data/simulated-epscs-1027.mat')
% load('/home/shababo/projects/mapping/code/psc-detection/data/simulated-data-longer-traces-epsc.mat')
clear timing_score
clear timing_score_tolerance
num_traces = length(event_times{1})

for i = 1:size(event_times,1)
    for j = 1:size(event_times,2)
        for k = 1:num_traces

        est_times = event_times{i,j}{k};
        [timing_score(i,j,k).correct_inds, timing_score(i,j,k).correct_err] = ...
            score_detection_trace(est_times, true_event_times{k}, tolerance);

        end
    end
end

timing_score_tolerance = tolerance;
save(resultsfile,'timing_score','timing_score_tolerance','-append')
    



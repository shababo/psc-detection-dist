function [map_curves, map_sample] = save_map(results_file)

load(results_file)
load(params.traces_filename)

if ~isfield(params,'start_ind')
    params.start_ind = 1;
end

if ~isfield(params,'duration')
    params.duration = size(traces,2);
end

traces = traces(:,params.start_ind:(params.start_ind + params.duration - 1));

if isfield(params,'traces_ind')
    traces = traces(params.traces_ind,:);
end

T = size(traces,2);

if isfield(params,'event_samples')
    event_samples = params.event_samples;
else
    event_samples = 4000;
end

map_curves = zeros(size(traces));
map_sample = struct();

for i = 1:size(traces,1)
    
    min_i = results(i).map_ind;
    this_curve = zeros(1,T);
    
    map_sample(i).tau = results(i).trials.tau{min_i};
    map_sample(i).time = results(i).trials.times{min_i};
    map_sample(i).amp = results(i).trials.amp{min_i};
    
    for j = 1:length(results(i).trials.times{min_i})
        
        
        ef = genEfilt_ar(results(i).trials.tau{min_i}{j},event_samples);
        [~, this_curve, ~] = addSpike_ar(results(i).trials.times{min_i}(j),...
                                            this_curve, 0, ef, ...
                                            results(i).trials.amp{min_i}(j),...
                                            results(i).trials.tau{min_i}{j},...
                                            traces(i,:), ...
                                            results(i).trials.times{min_i}(j), ...
                                            2, 1, 1);
                                        
    end
    map_curves(i,:) = this_curve + results(i).trials.base{min_i};
end

savename = [params.full_save_string(1:end-4) '-map.mat'];
save(savename,'map_curves','map_sample')
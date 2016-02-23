function map_trace = get_map_trace(results, params, T)

map_ind = results.map_ind;
    
map_trace = zeros(1,T);
       

for j = 1:length(results.trials.times{map_ind})

    ef = genEfilt_ar(results.trials.tau{map_ind}{j},params.event_samples);
    [~, map_trace, ~] = addSpike_ar(results.trials.times{map_ind}(j),...
                                        map_trace, 0, ef, ...
                                        results.trials.amp{map_ind}(j),...
                                        results.trials.tau{map_ind}{j},...
                                        zeros(1,T), ...
                                        results.trials.times{map_ind}(j), ...
                                        2, 1, 1);

end


map_trace = params.event_sign*(map_trace);% + results.trials.base{map_ind});
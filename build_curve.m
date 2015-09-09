function curve = build_curve(results,trace_ind, num_samples)

if isfield(results(1).params,'fBins')
    fBins = results(1).params.fBins;
else
    fBins = 4000;
end

min_i = results(trace_ind).min_err_ind;
    
    this_curve = zeros(1,num_samples);
    
    for j = 1:length(results(trace_ind).trials.times{min_i})
        
        ef = genEfilt_ar(results(trace_ind).trials.tau{min_i}{j},fBins);
        [~, this_curve, ~] = addSpike_ar(results(trace_ind).trials.times{min_i}(j),...
                                            this_curve, 0, ef, ...
                                            results(trace_ind).trials.amp{min_i}(j),...
                                            results(trace_ind).trials.tau{min_i}{j},...
                                            zeros(size(this_curve)), ...
                                            results(trace_ind).trials.times{min_i}(j), ...
                                            2, 1, 1);
                                        
    end
curve = -this_curve + results(trace_ind).trials.base{min_i};
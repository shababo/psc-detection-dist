function map_curve = plot_curve(results, trace, curveInd, plotFormat)

    if isempty(curveInd)
        [results.min_err, results.min_err_ind] = min(results.trials.obj);
        min_i = results.min_err_ind;
        curveInd = min_i;
    end
    fBins = 2000;
    N = length(trace);
    this_curve = zeros(1,N);
    for j = 1:length(results.trials.times{curveInd})

        ef = genEfilt_ar(results.trials.tau{curveInd}{j},min(fBins,N-results.trials.times{curveInd}(j)));

        [~, this_curve, ~] = addSpike_ar(results.trials.times{curveInd}(j),...
                                            this_curve, 0, ef, ...
                                            results.trials.amp{curveInd}(j),...
                                            results.trials.tau{curveInd}{j},...
                                            trace, ...
                                            results.trials.times{curveInd}(j), ...
                                            2, 1, 1);

    end
    map_curve = this_curve + results.trials.base{curveInd};
    
    if nargin>3
        plot(map_curve,plotFormat)
    end
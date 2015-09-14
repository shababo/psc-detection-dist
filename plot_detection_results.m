function plot_detection_results(results,var_names)

traces = [];
event_times = [];
target_feature_mat = [];

max_tau1 = Inf;
max_tau2 = Inf;
min_amp = 0;
min_time = 0;
max_time = Inf;





i = 1;
for j = 1:length(results)
    
%     if ~any(j == [18 21 31 32])
        
        if isfield(results(j).trials,'curves')
            traces(i,:) = results(j).trials.curves{results(i).min_err_ind}(min_time:end);
        else
            traces(i,:) = build_curve(results,j, 17500);
        end
        
        taus = zeros(0,2);
        good_events = [];
        
        for k = 1:length(results(i).trials.tau{results(i).min_err_ind})
            if all(results(i).trials.tau{results(i).min_err_ind}{k} < [max_tau1 max_tau2]) && results(i).trials.times{results(i).min_err_ind}(k) > min_time && results(i).trials.times{results(i).min_err_ind}(k) < max_time && results(i).trials.amp{results(i).min_err_ind}(k) > min_amp
                good_events = [good_events k];
                taus = [taus; results(i).trials.tau{results(i).min_err_ind}{k}];
            end
        end
        
        
        times = results(i).trials.times{results(i).min_err_ind}(good_events);
        event_times = [event_times times];

        

        new_features = [results(i).trials.amp{results(i).min_err_ind}(good_events)' taus results(i).trials.times{results(i).min_err_ind}(good_events)'/20000];
        target_feature_mat = [target_feature_mat; new_features];
        
        i = i+1;
        
%     end
    
end


% figure

% plot_trace_stack(-traces,zeros(size(traces)),[],zeros(length(traces),3),[],size(traces,2)-1)

figure
hist(event_times/20000,15)

size(target_feature_mat)

figure
plotmatrix(target_feature_mat)





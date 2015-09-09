function [target_feature_mat, trace_times, all_times] = plot_detection_results(results,var_names)

traces = [];
event_times = [];
target_feature_mat = [];

max_tau1 = Inf;
max_tau2 = Inf;
min_amp = 5;
min_time = 1;
max_time = Inf;
% hot_time_min = .048;
% hot_time_max = .05;
hot_time_min = .0047;
hot_time_max = .007;
% hot_amp_min = 0;
hot_amp_min = 20;
hot_amp_max = Inf;


for j = 1:length(results)
    
%     if ~any(j == [18 21 31 32])
        
%         if isfield(results(j).trials,'curves')
%             traces(j,:) = results(j).trials.curves{results(j).min_err_ind}(min_time:end);
%         else
%             traces(j,:) = build_curve(results,j, .05*20000);
%         end
        
        taus = zeros(0,2);
        time_constants = zeros(0,2);
        good_events = [];
        
        for k = 1:length(results(j).trials.tau{results(j).min_err_ind})
            if ~isempty(results(j).trials.tau{results(j).min_err_ind}{k})
                if all(results(j).trials.tau{results(j).min_err_ind}{k} < [max_tau1 max_tau2]) && results(j).trials.times{results(j).min_err_ind}(k) > min_time && results(j).trials.times{results(j).min_err_ind}(k) < max_time && results(j).trials.amp{results(j).min_err_ind}(k) > min_amp
                    good_events = [good_events k];
                    taus = [taus; results(j).trials.tau{results(j).min_err_ind}{k}];
                    
                    tau = results(j).trials.tau{results(j).min_err_ind}{k};
                    this_curve = zeros(1,1000);
                    ef = genEfilt_ar(tau,800);
                    [~, this_curve, ~] = addSpike_ar(100,...
                        this_curve, 0, ef, ...
                        1,...
                        tau,...
                        zeros(size(this_curve)), ...
                        100, ...
                        2, 1, 1);
                    

                    rise_time = find(this_curve > .67*max(this_curve),1,'first') - 100;
                    decay_time = find(this_curve(rise_time:end) < .33*max(this_curve),1,'first') + rise_time - 100;
                    
                    time_constants = [time_constants; rise_time decay_time];
                    
                    
                end
            end
        end
        
        
        times = results(j).trials.times{results(j).min_err_ind}(good_events);
        event_times = [event_times times];
        

        new_features = [results(j).trials.amp{results(j).min_err_ind}(good_events)' time_constants results(j).trials.times{results(j).min_err_ind}(good_events)'/20000];
        target_feature_mat = [target_feature_mat; new_features];
        
        
%     end
    
end

groups = [];
for i = 1:size(target_feature_mat,1)
    
    this_amp = target_feature_mat(i,1);
    this_time = target_feature_mat(i,4);
    if this_time > hot_time_min && this_time < hot_time_max && ...
            this_amp > hot_amp_min && this_amp < hot_amp_max
        disp('hit')
        groups = [groups 2];
    else
        groups = [groups 1];
    end
end

% figure

% plot_trace_stack(-traces,zeros(size(traces)),[],zeros(length(traces),3),[],size(traces,2)-1)

% figure
% hist(event_times/20000,15)

figure
gplotmatrix(target_feature_mat,target_feature_mat,groups)

figure
plotmatrix(target_feature_mat,var_names)

%%


all_times = [];
trace_times = cell(length(results),1);
for i = 1:length(results)
    trace_times{i} = [];
    this_trial_results = results(i).trials.times;
    for j = 1:length(this_trial_results)
        trace_times{i} = [trace_times{i} this_trial_results{j}];
        all_times = [all_times this_trial_results{j}];
    end  
end

figure;
for i = 1:length(results)
    subplot(length(results),1,i)
    hist(trace_times{i},100:100:2100)
    
    axis off
end

figure;
hist(all_times,20)


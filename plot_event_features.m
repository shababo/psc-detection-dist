function [target_feature_mat, groups] = plot_event_features(results_file,plot_params)

load(results_file)

traces = [];
event_times = [];
target_feature_mat = [];

min_tau1 = plot_params.min_tau1;
min_tau2 = plot_params.min_tau2;
max_tau1 = plot_params.max_tau1;
max_tau2 = plot_params.max_tau2;



min_amp = plot_params.min_amp;
min_time = plot_params.min_time;
max_time = plot_params.max_time;
hot_time_min = plot_params.hot_time_min;
hot_time_max = plot_params.hot_time_max;
hot_amp_min = plot_params.hot_amp_min;
hot_amp_max = plot_params.hot_amp_max;
hot_min_tau2 = plot_params.hot_min_tau2;


for j = 1:length(results)

        taus = zeros(0,2);
        time_constants = zeros(0,2);
        good_events = [];
        
        for k = 1:length(results(j).trials.tau{results(j).map_ind})
            if ~isempty(results(j).trials.tau{results(j).map_ind}{k})
                if all(results(j).trials.tau{results(j).map_ind}{k} > [min_tau1 min_tau2]) && all(results(j).trials.tau{results(j).map_ind}{k} < [max_tau1 max_tau2]) && results(j).trials.times{results(j).map_ind}(k) > min_time && results(j).trials.times{results(j).map_ind}(k) < max_time && results(j).trials.amp{results(j).map_ind}(k) > min_amp
                    
                    good_events = [good_events k];
                    taus = [taus; results(j).trials.tau{results(j).map_ind}{k}];
                    
                    tau = results(j).trials.tau{results(j).map_ind}{k};
                    this_curve = zeros(1,800);
                    ef = genEfilt_ar(tau,400);
                    [~, this_curve, ~] = addSpike_ar(1000,...
                        this_curve, 0, ef, ...
                        1,...
                        tau,...
                        zeros(size(this_curve)), ...
                        1, ...
                        2, 1, 1);
                    

                    rise_time = find(this_curve > .67*max(this_curve),1,'first');
                    [~,this_max_time] = max(this_curve);
                    decay_time = find(this_curve(rise_time:end) < .33*max(this_curve),1,'first') + rise_time - this_max_time;
                    
                    time_constants = [time_constants; rise_time decay_time];% this_max_time];
                    
                    
                end
            end
        end
        
        
        times = results(j).trials.times{results(j).map_ind}(good_events);
        event_times = [event_times times];
        

        new_features = [results(j).trials.amp{results(j).map_ind}(good_events)' taus results(j).trials.times{results(j).map_ind}(good_events)'/20000];
        target_feature_mat = [target_feature_mat; new_features];

end

groups = [];

for i = 1:size(target_feature_mat,1)
    
    this_amp = target_feature_mat(i,1);
    this_time = target_feature_mat(i,end);
    if this_time > hot_time_min && this_time < hot_time_max && ...
            this_amp > hot_amp_min && this_amp < hot_amp_max && ...
            target_feature_mat(i,3) > hot_min_tau2
        groups = [groups 2];
    else
        groups = [groups 1];
    end
end

figure
gplotmatrix(target_feature_mat,target_feature_mat,groups)

figure
plotmatrix(target_feature_mat,params.feature_names)

%%
% 
% all_times = [];
% trace_times = cell(length(results),1);
% for i = 1:length(results)
%     trace_times{i} = [];
%     this_trial_results = results(i).trials.times;
%     for j = 1:length(this_trial_results)
%         trace_times{i} = [trace_times{i} this_trial_results{j}];
%         all_times = [all_times this_trial_results{j}];
%     end  
% end
% 
% figure;
% for i = 1:length(results)
%     subplot(length(results),1,i)
%     bar(hist(trace_times{i},100:100:2100) / results(1).params.nsweeps)
%     ylim([0 1])
%     axis off
% end
% 
% figure;
% hist(all_times,20)

%% scratch

% || ...
%             (this_time > .1519 && this_time < .1526) || ...
%             (this_time > .2521 && this_time < .2531)) 
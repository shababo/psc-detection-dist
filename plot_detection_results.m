function [target_feature_mat, trace_times, all_times] = plot_detection_results(results,var_names)

traces = [];
event_times = [];
target_feature_mat = [];

min_tau1 = 0;
min_tau2 = 0;
max_tau1 = Inf;
max_tau2 = Inf;
min_amp = 10;
min_time = 1;
max_time = Inf;%.015*20000;
hot_time_min = .05;
hot_time_max = .055;
% hot_time_min = .0047;
% hot_time_max = .007;
hot_amp_min = 0;
% hot_amp_min = 20;
hot_amp_max = Inf;
hot_min_tau2 = 0;%50;


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
                if all(results(j).trials.tau{results(j).min_err_ind}{k} > [min_tau1 min_tau2]) && all(results(j).trials.tau{results(j).min_err_ind}{k} < [max_tau1 max_tau2]) && results(j).trials.times{results(j).min_err_ind}(k) > min_time && results(j).trials.times{results(j).min_err_ind}(k) < max_time && results(j).trials.amp{results(j).min_err_ind}(k) > min_amp
                    good_events = [good_events k];
                    taus = [taus; results(j).trials.tau{results(j).min_err_ind}{k}];
                    
                    tau = results(j).trials.tau{results(j).min_err_ind}{k};
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
        
        
        times = results(j).trials.times{results(j).min_err_ind}(good_events);
        event_times = [event_times times];
        

        new_features = [results(j).trials.amp{results(j).min_err_ind}(good_events)' taus results(j).trials.times{results(j).min_err_ind}(good_events)'/20000];
        target_feature_mat = [target_feature_mat; new_features];
        
        
%     end
    
end
size(target_feature_mat)
groups = [];
for i = 1:size(target_feature_mat,1)
    
    this_amp = target_feature_mat(i,1);
    this_time = target_feature_mat(i,end);
    if this_time > hot_time_min && this_time < hot_time_max && ...
            this_amp > hot_amp_min && this_amp < hot_amp_max && ...
            target_feature_mat(i,3) > hot_min_tau2
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
    bar(hist(trace_times{i},100:100:2100) / results(1).params.nsweeps)
    ylim([0 1])
    axis off
end

figure;
hist(all_times,20)


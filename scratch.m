% function output = scratch(linescan)
zero_ind = 6000;
duration = 4000;
fs = 20000;

dirname = '~/projects/mapping/data/0810/';

files = {'0810_cell02.mat','0810_cell04.mat','0810_cell06.mat','0810_cell08.mat','0810_cellXX.mat'};
%,'0810_cell07.mat'
num_cells = length(files);

condition_labels = {'950','900','850','800','750'};

num_conditions = length(condition_labels);

condition_inds = {{[12:16 22:26 47:51],[17:21 43:46],[27:31],[32:36],[37:41]},...
              {[18:24],[25:29],[35:39],[40:44],[45:48]},...
              {[14:18],[19:23],[33:37],[38:42],[43:48]},...
              {[14:18],[19:23],[24:28],[29:34],[35:39]},...
              {[28:35],[36:41],[42:46],[47:52],[53:59]}};

                        %{[7:15],[16:20],[21:25],[26:30],[]},...
dont_plot = 0;
data_ch = 1;

% get cell averages
cell_means = cell(num_cells,num_conditions);
cell_stds = cell(num_cells,num_conditions);
cell_current_max = zeros(num_cells,num_conditions);
cell_charge = zeros(num_cells,num_conditions);
condition_means = cell(num_conditions,1);

for i = 1:num_cells
    
    this_cell_file = [dirname files{i}];
    
    for j = 1:num_conditions
        

        traces = get_sweeps(this_cell_file,data_ch,condition_inds{i}{j},dont_plot);
        cell_means{i,j} = mean(traces);
        cell_means{i,j} = cell_means{i,j} - cell_means{i,j}(zero_ind);
        cell_stds{i,j} = std(traces);
        cell_current_max(i,j) = max(-cell_means{i,j}(zero_ind:zero_ind + duration));
        cell_charge(i,j) = sum(-cell_means{i,j}(zero_ind:zero_ind + duration));

        if i == 1
            condition_means{j} = cell_means{i,j}/num_cells;
        else
            condition_means{j} = condition_means{j} + cell_means{i,j}/num_cells;
        end
        
    end
end


condition_stds = cell(num_condtions,1);
for j = 1:num_conditions
    temp = [];
    for i = 1:num_cells
        temp(i,:) = cell_means{i,j};
    end
    condition_stds{j} = std(temp);
end
%%
normalized_cell_current_max = bsxfun(@ldivide,max(cell_current_max,[],2),cell_current_max);
normalized_cell_charge = bsxfun(@ldivide,max(cell_charge,[],2),cell_charge);



%% plot condition means


time = (0:length(condition_means{1})-1)/fs - zero_ind/fs;

figure;
set(gcf,'defaultAxesColorOrder',jet(num_conditions))
colors = {'r','m','g','c','b'};
for j = num_conditions:-1:1
    shadedErrorBar(time,condition_means{j} - condition_means{j}(zero_ind),condition_stds{j}/sqrt(num_cells),colors{j},1)
    hold on
end
hold off
% legend(fliplr(condition_labels))

ylim([-140 20])
xlim([-.02 .1])

%% plot summary

figure
plot([950 900 850 800 750],mean(normalized_cell_current_max))

figure
plot([950 900 850 800 750],mean(normalized_cell_charge))

%%
[targe_feature_mat_spike3, hot_grouping_spike3] = plot_event_features('data/evoked-pscs-strong-results-0000.mat',get_plot_params);
%%
[targe_feature_mat_spike1, hot_grouping_spike1] = plot_event_features('data/evoked-pscs-strong-results-0000.mat',get_plot_params);
%%
[targe_feature_mat_spike2, hot_grouping_spike2] = plot_event_features('data/evoked-pscs-strong-results-0000.mat',get_plot_params);
%%
amp1 = targe_feature_mat_spike1(hot_grouping_spike1 == 2,1);
amp2 = targe_feature_mat_spike2(hot_grouping_spike2 == 2,1);
amp3 = targe_feature_mat_spike3(hot_grouping_spike3 == 2,1);

bin_centers = 0:5:60;
num_bins = 8;
figure;
% histogram(amp3,bin_centers);
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75)
% hold on
histogram([amp3; amp2],bin_centers);
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','b','EdgeColor','w','facealpha',0.75)
hold on
histogram(amp1,bin_centers);
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','g','EdgeColor','w','facealpha',0.75)

figure;
histogram([amp1; amp2; amp3],bin_centers)


figure;
histogram([targe_feature_mat_spike1(:,1)],bin_centers);

%%

figure;
histogram(targe_feature_mat_spike1(hot_grouping_spike1 == 1 & hot_grouping_spike2 == 1 & hot_grouping_spike3 == 1,1),bin_centers);
hold on
histogram([amp1; amp2; amp3],bin_centers)

%%

amp1 = targe_feature_mat(hot_grouping == 2,1);


bin_centers = 0:2.5:40;
num_bins = 15;
figure;
% histogram(amp3,bin_centers);
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75)
% hold on
histogram(amp1,num_bins);
%%

figure;
histogram([amp1; amp2; amp3],bin_centers)


figure;
histogram([targe_feature_mat_spike1(:,1); targe_feature_mat_spike2(:,1); targe_feature_mat_spike3(:,1)],bin_centers);



figure;
histogram([targe_feature_mat_spike1(hot_grouping_spike1 == 1,1); targe_feature_mat_spike2(hot_grouping_spike2 == 1,1); targe_feature_mat_spike3(hot_grouping_spike3 == 1,1)],bin_centers);
hold on
histogram([amp1; amp2; amp3],bin_centers)

%% save parameter files for cluster

p_spike = [1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2];
a_min = [0 1e-2 1e-1 1e0 1e1 1e2 5e-2 5e-1 5e0 5e1 5e2];
num_sweeps = [1e1 1e2 1e3 1e4];

savefile_basename = 'data/cluster-param-files/params-pspike-%0.0e-amin-%0.0e-num_sweeps-%0.0e.mat';

for i = 1:length(p_spike)
    for j = 1:length(a_min)
        for k = 1:length(num_sweeps)
            params.cluster = 1;
            params.p_spike = p_spike(i);
            params.a_min = a_min(j);
            params.num_sweeps = num_sweeps(k);
            params = get_params(params);
            savefile = sprintf(savefile_basename,params.p_spike,params.a_min,params.num_sweeps);
            save(strrep(savefile,'+',''),'params')
            clear params
        end
    end
end
    
%% plot ROC

figure
for i = 2:size(roc,1)
    plot(roc(i,:,end,1),roc(i,:,end,2),'-o'); hold on
end
hold off
legend
title('each line is diff a_{min}, within varies p_{spike}')
%%
figure
for i = 1:size(roc_deconv,2)
    plot(roc_deconv(:,i,1),roc_deconv(:,i,2),'-o'); hold on
%% plot ROC
end
hold off
title('deconv')
figure
for i = 1:size(roc_crit,2)
    plot(roc_crit(:,i,1),roc_crit(:,i,2),'-o'); hold on
end
hold off
title('crit')




%%
    figure
for i = 2:size(roc,1)
    plot(squeeze(roc(i,end,:,1)),squeeze(roc(i,end,:,2)),'-o'); hold on
end
hold off
legend
title('each line is diff p_spike, within varies a_min')

%% plot ROC

figure
for i = 1:size(roc,2)
    plot(roc(:,i,1),roc(:,i,2),'-o'); hold on
end
hold off
legend

%% plot noise examples noise model figure

% load('data/work/example_noise_traces_work.mat')

figure;
subplot(221)
plot_trace_stack(traces,20,zeros(3,3),'-',[.010 10])
title('Voltage Clamp Recordings')

subplot(222)
plot_trace_stack(ar2_sim_data,20,zeros(3,3),'-',[.010 10])
title('Simulated Noise From Fits - AR(2)')

subplot(223)
plot_trace_stack(ar6_sim_data,20,zeros(3,3),'-',[.010 10])
title('Simulated Noise From Fits - AR(6)')

subplot(224)
plot_trace_stack(ar10_sim_data,20,zeros(3,3),'-',[.010 10])
title('Simulated Noise From Fits - AR(10)')

%% plot individual traces from real data example - fig0
load('/home/shababo/projects/mapping/code/psc-detection/stimfit/real-fit-events.mat')
figure;
legend_names = cell(1,3);
% for i = 1:3
    t = 0:1:.015*20000;
    tau_decay = event1_fit_params.tau1*20; decay = exp(-t/tau_decay);
    tau_rise = event1_fit_params.tau2*20; rise = -exp(-t/tau_rise);
    % plot(decay); hold on; plot(rise)
    plot(-(decay + rise)/max(decay+rise)*14,'linewidth',2)
    hold on; 
    legend_names{1} = ['event 1'];
    
    tau_decay = event2_fit_params.tau1*20; decay = exp(-t/tau_decay);
    tau_rise = event2_fit_params.tau2*20; rise = -exp(-t/tau_rise);
    % plot(decay); hold on; plot(rise)
    plot(-(decay + rise)/max(decay+rise)*19,'linewidth',2)
    hold on; 
    legend_names{2} = ['event 2'];
    
    tau_decay = event3_fit_params.tau1*20; decay = exp(-t/tau_decay);
    tau_rise = event3_fit_params.tau2*20; rise = -exp(-t/tau_rise);
    % plot(decay); hold on; plot(rise)
    plot(-(decay + rise)/max(decay+rise)*16,'linewidth',2)
    hold on; 
    legend_names{3} = ['event 3'];
% end
hold on
legend(legend_names)
axis off

%%

filenames = {'12_1_slice1_cell1.mat',...
             '12_1_slice1_cell2.mat',...
             '12_1_slice2_cell1.mat',...
             '12_1_slice3_cell1.mat',...
             '12_1_slice3_cell2.mat',...
             '12_1_slice4_cell1.mat',...
             '12_1_slice4_cell2.mat',...
             '12_1_slice5_cell1.mat',...
             '12_1_slice5_cell2.mat',...
             '12_2_slice1_cell1.mat',...
             '12_2_slice2_cell1.mat',...
             '12_2_slice2_cell2.mat',...
             '12_3_slice1_cell1.mat',...
             '12_3_slice1_cell2.mat',...
             '12_3_slice1_cell3.mat',...
             '12_3_slice2_cell1.mat',...
             '12_3_slice2_cell2.mat',...
             '12_3_slice3_cell1.mat',...
             '12_3_slice4_cell1.mat',...
             '12_3_slice4_cell2.mat',...
             '12_3_slice5_cell1.mat',...
             '12_3_slice6_cell1.mat',...
             '12_3_slice6_cell2.mat'};

         
         %%
         
for i = 1:length(filenames)
    compute_relative_obj_position(filenames{i},[]);
end

         %%
clear all
% 
% filenames = {'12_1_slice1_cell1.mat',...
% '12_1_slice1_cell2.mat',...
% '12_1_slice3_cell1.mat',...
% '12_1_slice3_cell2.mat',...
% '12_1_slice4_cell2.mat',...
% '12_1_slice5_cell1.mat',...
% '12_1_slice5_cell2.mat',...
% '12_2_slice1_cell1.mat',...
% '12_3_slice1_cell2.mat'};
% 
% 
% run_count_id = [9
% 2
% 10
% 3
% 13
% 15
% 2
% 6
% 6];
% 
trial_ids = [-60 -45 -30 -15 0 15 30 45 60 0 0 0 0 0 0 0 0 0
             0 0 0 0 0 0 0 0 0  -60 -45 -30 -15 0 15 30 45 60
             0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
  run_count_id = 8;       
filenames = {'12_3_slice1_cell2.mat'};
% filenames = {...
% '12_3_slice2_cell1.mat',...
% '12_3_slice3_cell1.mat'};
% 
% run_count_id = [
% 13
% 7];
% 
% trial_ids = [
%   -64.0000         0         0
%   -56.0000         0         0
%   -48.0000         0         0
%   -40.0000         0         0
%   -32.0000         0         0
%   -24.0000         0         0
%   -16.0000         0         0
%    -8.0000         0         0
%    0 0 0
%     7.5000         0         0
%    15.0000         0         0
%    22.5000         0         0
%    30.0000         0         0
%    37.5000         0         0
%    45.0000         0         0
%    52.5000         0         0
%    60.0000         0         0
%    0 -64 0
%             0  -56.0000         0
%          0  -48.0000         0
%          0  -40.0000         0
%          0  -32.0000         0
%          0  -24.0000         0
%          0  -16.0000         0
%          0   -8.0000         0
%                   0         0         0
%          0    8.0000         0
%          0   16.0000         0
%          0   24.0000         0
%          0   32.0000         0
%          0   40.0000         0
%                   0   48.0000         0
%          0   56.0000         0
%          0   64.0000         0
%          ];


peak_currents_cells_by_trial = cell(length(filenames),size(trial_ids,1));
peak_currents_cells_by_trial_mean = zeros(length(filenames),size(trial_ids,1));
peak_currents_cells_by_trial_std = zeros(length(filenames),size(trial_ids,1));
peak_currents_trial_mean = zeros(size(trial_ids,1),1);
peak_currents_trial_std = zeros(size(trial_ids,1),1);

baseline_window = 20000*[.299 .300]; measure_window = 20000*[.301 .330];

for i = 1:length(filenames)
    
    load(['data/' filenames{i}])
    
    [traces, traces_metadata] = get_sweeps_dir('data',filenames{i},0,1,0,Inf,'run_count',run_count_id(i));
    traces = traces{1};
    
    params1.run_count = run_count_id(i);
    match_inds = match_trials(params1, traces_metadata{1});
    traces = traces(match_inds,:);
    temp = traces_metadata{1};
    traces_metadata = temp(match_inds);
    
    trial_types = zeros(size(traces,1),1);
    
    for j = 1:size(trial_ids,1)
        
        params2.relative_position = trial_ids(j,:);
        match_inds = match_trials(params2, traces_metadata);
        size(match_inds)
        if isempty(match_inds)
            ['data/' filenames{i}]
%             assignin('base','whatthe',traces_metadata)
        end
        trial_types(match_inds) = j;
        
        peak_currents_cells_by_trial{i,j} = get_current_amp(traces(match_inds,:),baseline_window,measure_window);
        peak_currents_cells_by_trial_mean(i,j) = mean(peak_currents_cells_by_trial{i,j});
        peak_currents_cells_by_trial_std(i,j) = std(peak_currents_cells_by_trial{i,j});
        
    end


end

for j = 1:length(trial_ids)
    
    peak_currents_trial_mean(j) = mean(peak_currents_cells_by_trial_mean(:,j));
    peak_currents_trial_std(j) = mean(peak_currents_cells_by_trial_std(:,j));
    
end

%%
positions = -60:15:60;
colors = lines(length(positions));
switch_ind = length(peak_currents_trial_std)/2;
figure; %h = plot(positions,peak_currents_cells_by_trial_mean(:,1:switch_ind)');
hold on;
for i = 1:length(filenames)
    for j = 1:switch_ind
        plot(positions,peak_currents_cells_by_trial_mean(i,1:switch_ind),'color',colors(i,:))
        scatter(positions(j)*ones(length(peak_currents_cells_by_trial{i,j}),1),peak_currents_cells_by_trial{i,j},[],repmat(colors(i,:),length(peak_currents_cells_by_trial{i,j}),1),'filled');
        hold on;
    end
end
title('x')
figure; %h = plot(positions,peak_currents_cells_by_trial_mean(:,switch_ind+1:end)');
hold on;
for i = 1:length(filenames)
    for j = switch_ind+1:size(trial_ids,1)
        plot(positions,peak_currents_cells_by_trial_mean(i,switch_ind+1:end),'color',colors(i,:))
        hold on
        scatter(positions(j-switch_ind)*ones(length(peak_currents_cells_by_trial{i,j}),1),peak_currents_cells_by_trial{i,j},[],repmat(colors(i,:),length(peak_currents_cells_by_trial{i,j}),1),'filled');
        hold on;
    end
end
title('y')


%% normalization

peak_currents_cells_by_trial_norm = cell(length(filenames),size(trial_ids,1));
peak_currents_cells_by_trial_mean_norm = zeros(length(filenames),size(trial_ids,1));
peak_currents_cells_by_trial_std_norm = zeros(length(filenames),size(trial_ids,1));
peak_currents_trial_mean_norm = zeros(size(trial_ids,1),1);
peak_currents_trial_std_norm = zeros(size(trial_ids,1),1);


for i = 1:length(filenames)
    
    for j = 1:size(trial_ids,1)

        peak_currents_cells_by_trial_norm{i,j} = ...
            peak_currents_cells_by_trial{i,j}/max(peak_currents_cells_by_trial_mean(i,:));
        peak_currents_cells_by_trial_mean_norm(i,j) = mean(peak_currents_cells_by_trial_norm{i,j});
        peak_currents_cells_by_trial_std_norm(i,j) = std(peak_currents_cells_by_trial_norm{i,j});
        
    end
end
        
for j = 1:length(trial_ids)
    
    peak_currents_trial_mean_norm(j) = mean(peak_currents_cells_by_trial_mean_norm(:,j));
    peak_currents_trial_std_norm(j) = mean(peak_currents_cells_by_trial_std_norm(:,j));
    
end

%%

figure;

switch_ind = length(peak_currents_trial_std)/2;

figure; errorbar(peak_currents_cells_by_trial_mean_norm(:,1:switch_ind)',peak_currents_cells_by_trial_std_norm(:,1:switch_ind)');
title('x')
figure; errorbar(peak_currents_cells_by_trial_mean_norm(:,switch_ind+1:end)',peak_currents_cells_by_trial_std_norm(:,switch_ind+1:end)');
title('y')

%%

figure; errorbar(positions,peak_currents_trial_mean(1:9)',peak_currents_trial_std(1:9)');
title('x')
figure; errorbar(positions,peak_currents_trial_mean(10:end)',peak_currents_trial_std(10:end)');
title('y')

figure; errorbar(positions,peak_currents_trial_mean_norm(1:9)',peak_currents_trial_std_norm(1:9)');
title('x norm')
figure; errorbar(positions,peak_currents_trial_mean_norm(10:end)',peak_currents_trial_std_norm(10:end)');
title('y norm')

%%

clear all12_3_slice1_cell2

load(['data/12_3_slice6_cell2.mat'])
    
[traces, traces_metadata] = get_sweeps_dir('data','12_3_slice6_cell2.mat',0,1,0,Inf,'run_count',8);
traces = traces{1};

params1.run_count = 8;
match_inds = match_trials(params1, traces_metadata{1});
traces = traces(match_inds,:);
temp = traces_metadata{1};
traces_metadata = temp(match_inds);

unique_values = get_unique_metadata_vals(traces_metadata,'relative_position');

%%
traces_by_location = cell(11,11);

start_ind = 20000*.280; end_ind = 20000*.380;

for i = 1:size(traces,1)
    
    ind1 = traces_metadata(i).relative_position(1)/10 + 6;
    ind2 = traces_metadata(i).relative_position(2)/10 + 6;

    traces_by_location{ind1,ind2} = [traces_by_location{ind1,ind2} traces(i,start_ind:end_ind)'];
    
end



%%
%%

trace_count = 0;
for i = 1:size(traces_by_location,1)
    for j = 1:size(traces_by_location,2)
        
        trace_sum = sum(traces_by_location{i,j}',1);
        
        trace_count = trace_count + size(traces_by_location{i,j}',1);
        
    end
end

traces_avg = trace_sum/trace_count;

figure; plot(traces_avg)
     



%%

    
           filenames = {'12_3_slice1_cell1.mat',...
             '12_3_slice1_cell2.mat',...
             '12_3_slice1_cell3.mat',...
             '12_3_slice2_cell1.mat',...
             '12_3_slice2_cell2.mat',...
             '12_3_slice3_cell1.mat',...
             '12_3_slice4_cell1.mat',...
             '12_3_slice4_cell2.mat',...
             '12_3_slice5_cell1.mat',...
             '12_3_slice6_cell1.mat',...
             '12_3_slice6_cell2.mat'};
    
    
    
    
    
    
    
    
    
    
    
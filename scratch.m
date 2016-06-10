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
for i = 1:size(roc,1)
    plot(roc(i,:,end,1),roc(i,:,end,2),'-o'); hold on
end
hold off
legend
title('each line is diff a_{min}, within varies p_{spike}')
%%
figure

for i = 1:size(roc,2)
    plot(roc(:,i,1),roc(:,i,2),'-o'); hold on
end
hold off
title('deconv')
% figure
% for i = 1:size(roc_crit,2)
%     plot(roc_crit(:,i,1),roc_crit(:,i,2),'-o'); hold on
% end
% hold off
% title('crit')





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
    plot(roc(3:end,i,1),roc(3:end,i,2),'-o'); hold on
end



hold off
legend
%% all rocs
figure
plot(roc_bayes(3:end,1,1),roc_bayes(3:end,1,2),'.-b','linewidth',2,'markersize',25); hold on

plot(roc_cb(3:end,1,1),roc_cb(3:end,1,2),'.-g','linewidth',2,'markersize',25); hold on
% plot(roc_deconv(:,1,1),roc_deconv(:,1,2),'.-m','linewidth',2,'markersize',25); hold on
plot(roc_wiener(3:end,1,1),roc_wiener(3:end,1,2),'.-r','linewidth',2,'markersize',25); hold on
plot(roc_rand(:,1,1),roc_rand(:,1,2),'--k','linewidth',1,'markersize',5);

hold off

%% false pos

figure
plot(threshold(3:end),roc_cb(3:end,1,1),'.-g','linewidth',2,'markersize',25)
hold on
% plot(threshold,roc_deconv(:,1,1),'.-m','linewidth',2,'markersize',25)
% hold on
plot(threshold(3:end),roc_wiener(3:end,1,1),'.-r','linewidth',2,'markersize',25)
hold on
plot(threshold(3:end),roc_bayes(3:end,1,1),'.-b','linewidth',2,'markersize',25)

%% true pos

figure
plot(threshold(3:end),roc_cb(3:end,1,2),'.-g','linewidth',2,'markersize',25)
hold on
% plot(threshold,roc_deconv(:,1,2),'.-m','linewidth',2,'markersize',25)
% hold on
plot(threshold(3:end),roc_wiener(3:end,1,2),'.-r','linewidth',2,'markersize',25)
hold on
plot(threshold(3:end),roc_bayes(3:end,1,2),'.-b','linewidth',2,'markersize',25)

%% ratio tp/fp


figure
plot(threshold,roc_cb(:,1,2)./(roc_cb(:,1,1) + 1),'.-g','linewidth',2,'markersize',25)
hold on
plot(threshold,roc_deconv(:,1,2)./(roc_deconv(:,1,1) + 1),'.-m','linewidth',2,'markersize',25)
hold on
plot(threshold,roc_wiener(:,1,2)./(roc_wiener(:,1,1) + 1),'.-r','linewidth',2,'markersize',25)
hold on
plot(threshold,roc_bayes(:,1,2)./(roc_bayes(:,1,1) + 1),'.-b','linewidth',2,'markersize',25)


%% ratio fp/tp


figure
plot(threshold,1./(roc_cb(:,1,2)./roc_cb(:,1,1) ),'.-g','linewidth',2,'markersize',25)
hold on
plot(threshold,1./(roc_deconv(:,1,2)./roc_deconv(:,1,1) ),'.-m','linewidth',2,'markersize',25)
hold on
plot(threshold,1./(roc_wiener(:,1,2)./roc_wiener(:,1,1) ),'.-r','linewidth',2,'markersize',25)
hold on
plot(threshold,1./(roc_bayes(:,1,2)./roc_bayes(:,1,1) ),'.-b','linewidth',2,'markersize',25)

%% temporal accuracy histogram

bayes_errs = [];
cb_errs = [];
wiener_errs = [];

for i = 1:10
    
    bayes_errs = [bayes_errs double(squeeze(timing_score_bayes(6,1,i).correct_err))];
    cb_errs = [cb_errs double(squeeze(timing_score_cb(6,1,i).correct_err))];
    wiener_errs = [wiener_errs double(squeeze(timing_score_wiener(7,1,i).correct_err))];
    
end

figure
subplot(131)
histogram(cb_errs,'FaceColor','g','Normalization','pdf')
hold on
plot(ones(2,1)*median(cb_errs),[0 .25],'--k'); hold off
title(['Template Matching, mean error: ' num2str(median(cb_errs)/20) 'msec'])
ylim([0 .25])
xlim([0 20])
set(gca,'xticklabel',{'0','.25','.50','.75','1.0'})
subplot(132)
histogram(wiener_errs,'FaceColor','r','Normalization','pdf')
hold on
plot(ones(2,1)*median(wiener_errs),[0 .25],'--k'); hold off
title(['Wiener Filter, mean error: ' num2str(median(wiener_errs)/20) 'msec'])
ylim([0 .25])
xlim([0 20])
set(gca,'xticklabel',{'0','.25','.50','.75','1.0'})
subplot(133)
histogram(bayes_errs,'FaceColor','b','Normalization','pdf')
hold on
plot(ones(2,1)*median(bayes_errs),[0 .25],'--k'); hold off
title(['Bayesian, mean error: ' num2str(median(bayes_errs)/20) 'msec'])
ylim([0 .25])
xlim([0 20])
set(gca,'xticklabel',{'0','.25','.50','.75','1.0'})

%% simulate random detection

rates = [0 .1 .5 1 2 4 8 16 32];

event_times = cell(length(rates),1);

for i = 1:length(rates)
    event_times{i} = cell(1,10);
    for j = 1:10
        num_events = poissrnd(rates(i));
        event_times{i}{j} = ceil(rand(1,num_events)*20000);
    end
end

save('data/random-detection-results-harder.mat','event_times');

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
edi
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
filenames = {'12_17_slice1chrims_cell1.mat',...
    '12_17_slice1chrims_cell2.mat',...
    '12_17_slice2chrims_cell1.mat',...
    '12_17_slice2chrims_cell2.mat',...
    '12_17_slice2chrims_cell3.mat',...
    '12_17_slice3chrims_cell1.mat',...
    '12_17_slice3chrims_cell2.mat',...
    '12_17_slice4chrims_cell1.mat',...
    '12_17_slice4chrims_cell2.mat',...
    '12_17_slice5chrims_cell1.mat',...
    '12_17_slice5chrims_cell2.mat'};
         
         %%
         
for i = 1:length(filenames)
    compute_relative_obj_position(['data/' filenames{i}],[]);
end

%% Chrimson Good Currents
clear all
% 
filenames = {...%'12_1_slice1_cell1.mat',...
'12_1_slice1_cell2.mat',...
'12_1_slice3_cell1.mat',...
...%'12_1_slice3_cell2.mat',...
'12_1_slice4_cell2.mat',...
'12_1_slice5_cell1.mat',...
'12_1_slice5_cell2.mat',...
'12_2_slice1_cell1.mat',...
'12_3_slice1_cell2.mat',...
...%'12_3_slice2_cell1.mat',...
...% '12_3_slice3_cell1.mat',...
'12_17_slice1chrims_cell1.mat',...
'12_17_slice1chrims_cell2.mat',...
'12_17_slice2chrims_cell1.mat',...
'12_17_slice3chrims_cell1.mat',...
'12_17_slice3chrims_cell2.mat'};


run_count_id = {...%9
2,...
10,...
...%3
13,...
15,...
2,...
6,...
6,...
...%13,...
...%7,...
{3,4},{5},{4,5},{5,6},{4,5}};

%% Chrimson w/ TF

filenames = {'12_19_slice1_cell1.mat',...
    '12_19_slice1_cell2.mat',...
    '12_19_slice1_cell3.mat',...
    '12_20_slice2_cell3.mat'};

run_count_id = {3, 3, 3 ...
    3};

trial_ids1 = [-60 -45 -30 -15 0 15 30 45 60 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
             0 0 0 0 0 0 0 0 0  -60 -45 -30 -15 0 15 30 45 60 0 0 0 0 0 0 0 0 0
             0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -60 -45 -30 -15 0 15 30 45 60]';
         
%% Chrimson w/o TF XYZ

filenames = {'12_22_slice2_cell1.mat',...
    '12_22_slice4_cell1.mat',...
    '12_22_slice4_cell2.mat',...
    '12_22_slice4_cell4.mat'};

run_count_id = {3, {3,4}, 7 ...
    2};

trial_ids1 = [-60 -45 -30 -15 0 15 30 45 60 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
             0 0 0 0 0 0 0 0 0  -60 -45 -30 -15 0 15 30 45 60 0 0 0 0 0 0 0 0 0
             0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -60 -45 -30 -15 0 15 30 45 60]';
         
%% soma-chr2 good currents
clear all
filenames = {'12_15_slice1_cell2.mat','12_15_slice1_cell3.mat','12_15_slice2_cell1.mat','12_15_slice3_cell1.mat',...
    '12_17_slice1_cell1.mat','12_17_slice1_cell3.mat','12_17_slice2_cell1.mat','12_17_slice3_cell3.mat','12_17_slice3_cell4.mat',}

  run_count_id = {4, 3, {4, 5}, {3, 4},{5,6,7},{3,4},{3,4,5},{3,4},4}; 



  
%   run_count_id = 8;
% filenames = {'12_3_slice1_cell2.mat'};
% filenames = {...
% '12_3_slice2_cell1.mat',...
% '12_3_slice3_cell1.mat'};
% 
% run_count_id = [
% 13
% 7];
% 

%% trial ids

trial_ids1 = [-60 -45 -30 -15 0 15 30 45 60 0 0 0 0 0 0 0 0 0 
             0 0 0 0 0 0 0 0 0  -60 -45 -30 -15 0 15 30 45 60
             0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';

trial_ids2 = [
  -64.0000         0         0
  -48.0000         0         0
  -32.0000         0         0
  -16.0000         0         0
   0 0 0
   15.0000         0         0
   30.0000         0         0
   45.0000         0         0
   60.0000         0         0
   0 -64 0
         0  -48.0000         0
         0  -32.0000         0
         0  -16.0000         0
                  0         0         0
         0   16.0000         0
         0   32.0000         0
                  0   48.0000         0
         0   64.0000         0
         ];
     
trial_ids3 = [
  -64.0000         0         0
  -56.0000         0         0
  -48.0000         0         0
  -40.0000         0         0
  -32.0000         0         0
  -24.0000         0         0
  -16.0000         0         0
   -8.0000         0         0
   0 0 0
    7.5000         0         0
   15.0000         0         0
   22.5000         0         0
   30.0000         0         0
   37.5000         0         0
   45.0000         0         0
   52.5000         0         0
   60.0000         0         0
   0 -64 0
            0  -56.0000         0
         0  -48.0000         0
         0  -40.0000         0
         0  -32.0000         0
         0  -24.0000         0
         0  -16.0000         0
         0   -8.0000         0
                  0         0         0
         0    8.0000         0
         0   16.0000         0
         0   24.0000         0
         0   32.0000         0
         0   40.0000         0
                  0   48.0000         0
         0   56.0000         0
         0   64.0000         0
         ];

%%

peak_currents_cells_by_trial = cell(length(filenames),size(trial_ids1,1));
peak_currents_cells_by_trial_mean = zeros(length(filenames),size(trial_ids1,1));
peak_currents_cells_by_trial_std = zeros(length(filenames),size(trial_ids1,1));
peak_currents_trial_mean = zeros(size(trial_ids1,1),1);
peak_currents_trial_std = zeros(size(trial_ids1,1),1);

baseline_window = 20000*[.299 .300]; measure_window = 20000*[.301 .330];

for i = 1:length(filenames)
    
    load(['data/' filenames{i}])
    
    [traces, traces_metadata] = get_sweeps_dir('data',filenames{i},0,1,0,Inf,'run_count',run_count_id{i});
    traces = traces{1};
    
    params1.run_count = run_count_id{i};
    match_inds = match_trials(params1, traces_metadata{1});
    traces = traces(match_inds,:);
    temp = traces_metadata{1};
    traces_metadata = temp(match_inds);
    
    trial_types = zeros(size(traces,1),1);
    

    for j = 1:size(trial_ids1,1)
        
        params2.relative_position = trial_ids1(j,:);
%         params3.relative_position = trial_ids2(j,:);
%         match_inds = unique([match_trials(params2, traces_metadata) match_trials(params3, traces_metadata)]);
        match_inds = unique([match_trials(params2, traces_metadata)]);
        size(match_inds)
        if isempty(match_inds)
            ['data/' filenames{i}]
%             assignin('base','whatthe',traces_metadata)
        end
        trial_types(match_inds) = j;
        
        upper_limit = 500;
        
        peak_currents_cells_by_trial{i,j} = get_current_amp(traces(match_inds,:),baseline_window,measure_window);
        peak_currents_cells_by_trial{i,j}(peak_currents_cells_by_trial{i,j} > upper_limit) = [];
        peak_currents_cells_by_trial_mean(i,j) = mean(peak_currents_cells_by_trial{i,j});
        peak_currents_cells_by_trial_std(i,j) = std(peak_currents_cells_by_trial{i,j});
        
    end


end


for j = 1:length(trial_ids1)
    
    peak_currents_trial_mean(j) = mean(peak_currents_cells_by_trial_mean(:,j));
    peak_currents_trial_std(j) = mean(peak_currents_cells_by_trial_std(:,j));
    
end

%% X Y ONLY
positions = -60:15:60;
colors = lines(length(filenames));
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
% legend(filenames)
title('x')
figure; %h = plot(positions,peak_currents_cells_by_trial_mean(:,switch_ind+1:end)');
hold on;
for i = 1:length(filenames)
    for j = switch_ind+1:size(trial_ids1,1)
        plot(positions,peak_currents_cells_by_trial_mean(i,switch_ind+1:end),'color',colors(i,:))
        hold on
        scatter(positions(j-switch_ind)*ones(length(peak_currents_cells_by_trial{i,j}),1),peak_currents_cells_by_trial{i,j},[],repmat(colors(i,:),length(peak_currents_cells_by_trial{i,j}),1),'filled');
        hold on;
    end
end
% legend(filenames)
title('y')
%%
positions = -60:15:60;

%% X Y Z


colors = lines(length(filenames));
switch_ind = length(peak_currents_trial_std)/3;
figure; %h = plot(positions,peak_currents_cells_by_trial_mean(:,1:switch_ind)');
hold on;
for i = 1:length(filenames)
    for j = 1:switch_ind
        plot(positions,peak_currents_cells_by_trial_mean(i,1:switch_ind),'color',colors(i,:))
        scatter(positions(j)*ones(length(peak_currents_cells_by_trial{i,j}),1),peak_currents_cells_by_trial{i,j},[],repmat(colors(i,:),length(peak_currents_cells_by_trial{i,j}),1),'filled');
        hold on;
    end

end
% legend(filenames)
title('x')
figure; %h = plot(positions,peak_currents_cells_by_trial_mean(:,switch_ind+1:end)');
hold on;
for i = 1:length(filenames)
    for j = switch_ind+1:switch_ind*2
        plot(positions,peak_currents_cells_by_trial_mean(i,switch_ind+1:switch_ind*2),'color',colors(i,:))
        hold on
        scatter(positions(j-switch_ind)*ones(length(peak_currents_cells_by_trial{i,j}),1),peak_currents_cells_by_trial{i,j},[],repmat(colors(i,:),length(peak_currents_cells_by_trial{i,j}),1),'filled');
        hold on;
    end
end
% legend(filenames)
title('y')

figure; %h = plot(positions,peak_currents_cells_by_trial_mean(:,switch_ind+1:end)');
hold on;
for i = 1:length(filenames)
    for j = switch_ind*2+1:size(trial_ids1,1)
        plot(positions,peak_currents_cells_by_trial_mean(i,switch_ind*2+1:end),'color',colors(i,:))
        hold on
        scatter(positions(j-switch_ind*2)*ones(length(peak_currents_cells_by_trial{i,j}),1),peak_currents_cells_by_trial{i,j},[],repmat(colors(i,:),length(peak_currents_cells_by_trial{i,j}),1),'filled');
        hold on;
    end
end
% legend(filenames)
title('z')
%% normalization

peak_currents_cells_by_trial_norm = cell(length(filenames),size(trial_ids1,1));
peak_currents_cells_by_trial_mean_norm = zeros(length(filenames),size(trial_ids1,1));
peak_currents_cells_by_trial_std_norm = zeros(length(filenames),size(trial_ids1,1));
peak_currents_trial_mean_norm = zeros(size(trial_ids1,1),1);
peak_currents_trial_std_norm = zeros(size(trial_ids1,1),1);


for i = 1:length(filenames)
    
    for j = 1:size(trial_ids1,1)
        
        start = floor((j-1)/9)*9 + 1
        stop = start + 8
        peak_currents_cells_by_trial_norm{i,j} = ...
            peak_currents_cells_by_trial{i,j}/max(peak_currents_cells_by_trial_mean(i,start:stop));
        peak_currents_cells_by_trial_mean_norm(i,j) = mean(peak_currents_cells_by_trial_norm{i,j});
        peak_currents_cells_by_trial_std_norm(i,j) = std(peak_currents_cells_by_trial_norm{i,j});
        
    end
end
        
for j = 1:length(trial_ids1)
    
    peak_currents_trial_mean_norm(j) = mean(peak_currents_cells_by_trial_mean_norm(:,j));
    peak_currents_trial_std_norm(j) = mean(peak_currents_cells_by_trial_std_norm(:,j));
    
end

%%


switch_ind = length(peak_currents_trial_std)/2;

figure; errorbar(peak_currents_cells_by_trial_mean_norm(:,1:switch_ind)',peak_currents_cells_by_trial_std_norm(:,1:switch_ind)');
title('x')
figure; errorbar(peak_currents_cells_by_trial_mean_norm(:,switch_ind+1:end)',peak_currents_cells_by_trial_std_norm(:,switch_ind+1:end)');
title('y')

%%


switch_ind = length(peak_currents_trial_std)/3;

figure; errorbar(peak_currents_cells_by_trial_mean_norm(:,1:switch_ind)',peak_currents_cells_by_trial_std_norm(:,1:switch_ind)');
title('x')
figure; errorbar(peak_currents_cells_by_trial_mean_norm(:,switch_ind+1:switch_ind*2)',peak_currents_cells_by_trial_std_norm(:,switch_ind+1:switch_ind*2)');
title('y')
figure; errorbar(peak_currents_cells_by_trial_mean_norm(:,switch_ind*2+1:end)',peak_currents_cells_by_trial_std_norm(:,switch_ind*2+1:end)');
title('z')

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

figure; errorbar(positions,peak_currents_trial_mean(1:9)',peak_currents_trial_std(1:9)');
title('x')
figure; errorbar(positions,peak_currents_trial_mean(10:end)',peak_currents_trial_std(10:end)');
title('y')

%%
figure; errorbar(positions,peak_currents_trial_mean_norm(1:9)',peak_currents_trial_std_norm(1:9)');
title('x norm')
figure; errorbar(positions,peak_currents_trial_mean_norm(10:18)',peak_currents_trial_std_norm(10:18)');
title('y norm')
figure; errorbar(positions,peak_currents_trial_mean_norm(19:end)',peak_currents_trial_std_norm(19:end)');
title('z norm')

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
         
         
         
%% testing wiener filter

trace = traces{1,1}(1,:);
nfft = length(trace) + length(template);
dt = 1/20000;
threshold = 2.0;
min_window = 20;

ar_noise_params.sigma_sq = 2.948727926352792;
ar_noise_params.phi = [1.000000000000000, -0.982949319747574, 0.407063852831604];
gamma = 1e6;

figure;
subplot(1,2,1)
[filtered_trace, event_times_init] = wiener_filter(trace,template,ar_noise_params,...
            nfft, 1/20000, threshold, min_window);
        
subplot(1,2,2)
plot(filtered_trace/max(filtered_trace)*5)
hold on
plot(trace + 40 - trace(1))
hold on
scatter(event_times_init, 8*ones(1,length(event_times_init)))

%% spike detection on grid
 
trace_grids_3_31_s2c2_r2_3 = {traces_by_location_3_31_s2c2_r3_5mw, traces_by_location_3_31_s2c2_r3_10mw, traces_by_location_3_31_s2c2_r3_15mw,...
    traces_by_location_3_31_s2c2_r3_25mw, traces_by_location_3_31_s2c2_r3_50mw, traces_by_location_3_31_s2c2_r3_100mw};

trace_grids_3_29_s1c2_r2 = {traces_by_location_3_29_s1c2_r2_25mw, traces_by_location_3_29_s1c2_r2_50mw, traces_by_location_3_29_s1c2_r2_100mw};

trace_grids_3_29_s1c4_r2 = {traces_by_location_3_29_s1c4_r2_25mw, traces_by_location_3_29_s1c4_r2_50mw, traces_by_location_3_29_s1c4_r2_100mw};

trace_grids_3_31_s1c1_r4_5 = {traces_by_location_3_31_s1c1_r4_5_25mw, traces_by_location_3_31_s1c1_r4_5_50mw, traces_by_location_3_31_s1c1_r4_5_100mw};

trace_grids_3_31_s1c2_r4_5 = {traces_by_location_3_31_s1c2_r5_5mw, traces_by_location_3_31_s1c2_r5_10mw, traces_by_location_3_31_s1c2_r5_15mw,...
    traces_by_location_3_31_s1c2_r4_25mw, traces_by_location_3_31_s1c2_r4_50mw, traces_by_location_3_31_s1c2_r4_100mw};

trace_grids_4_5_s2c1_r5 = {traces_by_location_4_5_s2c1_r5_25mw, traces_by_location_4_5_s2c1_r5_50mw, traces_by_location_4_5_s2c1_r5_100mw};

trace_grids_4_6_s3c2_r1 = {traces_by_location_4_6_s3c2_r1_25mw, traces_by_location_4_6_s3c2_r1_50mw, traces_by_location_4_6_s3c2_r1_100mw};

traces_trids_4_6_s3c5_r1 = {traces_by_location_4_6_s3c5_r1_25mw, traces_by_location_4_6_s3c5_r1_50mw, traces_by_location_4_6_s3c5_r1_100mw};

traces_trids_4_6_s3c7_r2 = {traces_by_location_4_6_s3c7_r2_25mw, traces_by_location_4_6_s3c7_r2_50mw, traces_by_location_4_6_s3c7_r2_100mw};

traces_trids_4_6_s3c8_r3 = {traces_by_location_4_6_s3c8_r3_25mw, traces_by_location_4_6_s3c8_r3_50mw, traces_by_location_4_6_s3c8_r3_100mw};

%%

trace_grids = traces_trids_4_6_s3c8_r3;

detection_grids = cell(size(trace_grids));

for i = 1:length(trace_grids)
    
    trace_grid_tmp = trace_grids{i};
    [traces_tmp, rebuild_map] = stack_traces(trace_grid_tmp);

    detection_results = detect_peaks(-1.0*bsxfun(@minus,traces_tmp,median(traces_tmp,2)),4.0,20,1,1,0)*70;
    detection_grids{i} = unstack_traces(detection_results,rebuild_map);
    
end

detection_results_4_6_s3c8_r3 = detection_results;
detection_grids_4_6_s3c8_r3 = detection_grids;

%%


figure; compare_trace_stack_grid({trace_grids{:},detection_grids_4_6_s3c8_r3{:}},...
    5,1,0,{'25 mW', '50 mW', '100 mW'},2)


%% count spikes and get means

spike_counts = zeros([size(detection_grids{1}) length(detection_grids)]);
max_val = 0;

axs = [];

figure
colormap hot
for i = 1:length(detection_grids)
    
    this_grid = detection_grids{i};
    for j = 1:size(this_grid,1)
        for k = 1:size(this_grid,2)
            
            spike_counts(j,k,i) = length(find(this_grid{j,k}(:)))/size(this_grid{j,k},1);
            if spike_counts(j,k,i) > max_val
                max_val = spike_counts(j,k,i);
            end
            
        end
    end
    subplot(1,length(detection_grids),i)
%     subplot(2,ceil(length(detection_grids)/2),i)
    imagesc(spike_counts(:,:,i));
    axis square
    axis off
end

for i = 1:length(detection_grids)
    subplot(1,length(detection_grids),i)
%     subplot(2,ceil(length(detection_grids)/2),i)
    caxis([0 max_val])
end

spike_counts_4_6_s3c8_r3 = spike_counts;

%% delay times and get means

delays = zeros([size(detection_grids{1}) length(detection_grids)]);
max_val = 0;

axs = [];

figure

colormap hot
for i = 1:length(detection_grids)
    
    this_grid = detection_grids{i};
    for j = 1:size(this_grid,1)
        for k = 1:size(this_grid,2)
            
            these_delays = [];
            
            for m = 1:size(this_grid{j,k},1)
                these_delays = [these_delays find(this_grid{j,k}(m,:),1,'first')/20000 - .005];
            end
            
            delays(j,k,i) = mean(these_delays);
            
            if delays(j,k,i) > max_val
                max_val = delays(j,k,i);
            end
            
        end
    end
    subplot(2,ceil(length(detection_grids)/2),i)
    pcolor(delays(:,:,i));
    axis ij
    axis square
    axis off
end

delays_4_6_s3c8_r3 = delays;

for i = 1:length(detection_grids)
    subplot(2,ceil(length(detection_grids)/2),i)
    caxis([0 max_val])
end

%%

figure 

for i = 1:3
    
    subplot(3,3,i)
    pcolor(flipud(delays_3_31_s2c2_r2_3(:,:,i+3)));
    caxis([0 .05])
    axis square
    axis off
end

for i = 1:3
    
    subplot(3,3,i+3)
    pcolor(flipud(delays_3_29_s1c2_r2(:,:,i)));
    caxis([0 .05])
    axis square
    axis off
end


for i = 1:3
    
    subplot(3,3,i+6)
    pcolor(flipud(delays_3_31_s1c2_r4_5(:,:,i+3)));
    caxis([0 .05])
    axis square
    axis off
end


colormap hot

%%

figure 

for i = 1:3
    
    subplot(3,3,i)
    imagesc(spike_counts_4_6_s3c5_r1(:,:,i));
    caxis([0 2])
    axis square
    axis off
end

for i = 1:3
    
    subplot(3,3,i+3)
    imagesc(spike_counts_4_6_s3c2_r1(:,:,i));
    caxis([0 2])
    axis square
    axis off
end


for i = 1:3
    
    subplot(3,3,i+6)
    imagesc(spike_counts_4_6_s3c7_r2(:,:,i));
    caxis([0 2])
    axis square
    axis off
end
colormap hot

%%
figure
for i = 1:6
    
    subplot(2,6,i)
    imagesc(spike_counts_3_31_s2c2_r2_3(:,:,i));
    caxis([0 3])
%     axis square
    axis off
end

for i = 1:6
    
    subplot(2,6,i + 6)
    imagesc(spike_counts_3_31_s1c2_r4_5(:,:,i));
    caxis([0 3])
%     axis square
    axis off
end



colormap hot

%%

%%
trace_grid_ch1 = traces_by_location_5_12_s2c1_2_r4{1};
trace_grid_ch2= traces_by_location_5_12_s2c1_2_r4{2};

across_ch_corr_image = zeros(size(trace_grid_ch1));

for i = 1:size(trace_grid_ch1,1)
    for j = 1:size(trace_grid_ch1,2)
        for k = 1:size(trace_grid_ch1{i,j},1)
            corr_mat = corr([trace_grid_ch1{i,j}(k,:); trace_grid_ch2{i,j}(k,:)]');
            across_ch_corr_image(i,j) = across_ch_corr_image(i,j) + corr_mat(1,2);
        end
        across_ch_corr_image(i,j) = across_ch_corr_image(i,j)/size(trace_grid_ch1{i,j},1);
    end
end

figure; 
imagesc(across_ch_corr_image)
colormap hot
colorbar

%%

event_timeseries2 = get_event_times_init(results2,2000,1,10);
event_timeseries1 = get_event_times_init(results,2000,1,10);
event_timeseries1_smooth = smoothts(event_timeseries1,'g',100,20);
event_timeseries2_smooth = smoothts(event_timeseries2,'g',100,20);
events_ts_grid1_smooth = unstack_traces(event_timeseries1_smooth*1000,params.rebuild_map);
events_ts_grid2_smooth = unstack_traces(event_timeseries2_smooth*1000,params2.rebuild_map);
figure; compare_trace_stack_grid_overlap({events_ts_grid1_smooth,events_ts_grid2_smooth},3,1,[],0,{'L4','L5'},1)

%%

max_xcorr = cell(size(events_ts_grid1_smooth));
mad_xcorr_lag = cell(size(events_ts_grid1_smooth));

for i = 1:size(events_ts_grid1_smooth,1)
    for j = 1:size(events_ts_grid1_smooth,2)
        max_xcorr{i,j} = zeros(size(events_ts_grid1_smooth{i,j},1),1);
        mad_xcorr_lag{i,j} = zeros(size(events_ts_grid1_smooth{i,j},1),1);
        for k = 1:size(events_ts_grid1_smooth{i,j},1)
            [max_xcorr{i,j}(k), mad_xcorr_lag{i,j}(k)] = max(xcorr(events_ts_grid1_smooth{i,j}(k,:),events_ts_grid1_smooth{i,j}(k,:)));
        end
    end
end

%%

max_xcorr_img = zeros(size(events_ts_grid1_smooth));
mad_xcorr_lag_img = zeros(size(events_ts_grid1_smooth));

for i = 1:size(events_ts_grid1_smooth,1)
    for j = 1:size(events_ts_grid1_smooth,2)
        max_xcorr_img(i,j) = mean(max_xcorr{i,j});
        mad_xcorr_lag_img(i,j) = mean(mad_xcorr_lag{i,j});

    end
end

figure;
subplot(121)
imagesc(max_xcorr_img)
colorbar
subplot(122)
imagesc(mad_xcorr_lag_img)
colorbar

%%


xcorrs = cell(size(events_ts_grid1_smooth));

for i = 1:size(events_ts_grid1_smooth,1)
    for j = 1:size(events_ts_grid1_smooth,2)
        xcorrs{i,j} = zeros(size(events_ts_grid1_smooth{i,j}));
        for k = 1:size(events_ts_grid1_smooth{i,j},1)
            xcorrs{i,j}(k,:) = max(xcorr(events_ts_grid1_smooth{i,j}(k,:),events_ts_grid1_smooth{i,j}(k,:)));
        end
    end
end






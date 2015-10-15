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
for i = 1:size(roc,2)
    plot(roc(2:end,i,end,1),roc(2:end,i,end,2),'-o'); hold on
end
hold off
legend
title('each line is diff p_spike, within varies a_min')


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

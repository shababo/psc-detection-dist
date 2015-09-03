% function output = scratch(linescan)
zero_ind = 6000;
duration = 400;
fs = 20000;

dirname = '~/Desktop/0810/';

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

%%
condition_stds = cell(num_conditions,1);
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
%     shadedErrorBar(time,condition_means{j} - condition_means{j}(zero_ind),condition_stds{j}/sqrt(num_cells),colors{j},1)
    plot(time,condition_means{j} - condition_means{j}(zero_ind),colors{j},'linewidth',2)
    hold on
end
hold off
legend(fliplr(condition_labels))

xlabel('time (sec)')
ylabel('current (pA)')


ylim([-140 20])
xlim([-.02 .1])

%% plot summary

wavelengths = [950 900 850 800 750];
figure
for i = 1:num_cells
    plot(wavelengths,normalized_cell_current_max(i,:),'bo-');
    hold on;
end
scatter([950 900 850 800 750],mean(normalized_cell_current_max),75,'b','filled')
%%
figure
for i = 1:length(wavelengths)
    scatter(wavelengths(i)*ones(size(normalized_cell_charge,1),1),normalized_cell_charge(:,i),40,'b');
    hold on;
end
scatter([950 900 850 800 750],mean(normalized_cell_charge),75,'b','filled')
%%


figure
for i = 1:num_cells
    plot(wavelengths,cell_current_max(i,:),'bo-');
    hold on;
end
scatter([950 900 850 800 750],mean(cell_current_max),75,'b','filled')


    
    
    
    
    
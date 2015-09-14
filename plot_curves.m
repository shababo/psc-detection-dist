function plot_curves(traces_file, results_file)

load(traces_file)

load(results_file)
traces = traces(:,results(1).params.start_ind:(results(1).params.start_ind + results(1).params.duration));



T = size(traces,2);
fBins = 2000;

map_curves = zeros(size(traces));

offset = 50;


for i = 1:size(traces,1)
    
    min_i = results(i).min_err_ind;
    this_curve = zeros(1,T);
    
    for j = 1:length(results(i).trials.times{min_i})
        
        ef = genEfilt_ar(results(i).trials.tau{min_i}{j},fBins);
        [~, this_curve, ~] = addSpike_ar(results(i).trials.times{min_i}(j),...
                                            this_curve, 0, ef, ...
                                            results(i).trials.amp{min_i}(j),...
                                            results(i).trials.tau{min_i}{j},...
                                            traces(i,:), ...
                                            results(i).trials.times{min_i}(j), ...
                                            2, 1, 1);
                                        
    end
    map_curves(i,:) = -this_curve + results(i).trials.base{min_i};
end


figure;
plot_trace_stack(traces,zeros(size(traces)),[],bsxfun(@plus,zeros(length(traces),3),[1 .4 .4]),[],size(traces,2)-1,offset)
hold on
plot_trace_stack(map_curves,zeros(size(traces)),[],bsxfun(@plus,zeros(length(traces),3),[0 0 1]),[],size(traces,2)-1,offset)

% % figure;
% % traces = traces(2:4:end,:);
% plot_trace_stack(traces,zeros(size(traces)),[],bsxfun(@plus,zeros(length(traces),3),[.4 .4 1]),[],size(traces,2)-1,50)
% hold on
% % map_curves = map_curves(2:4:end,:);
% plot_trace_stack(map_curves,zeros(size(traces)),[],bsxfun(@plus,zeros(length(traces),3),[0 0 1]),[],size(traces,2)-1,50)
% 
% title_string = ['a-min = ' num2str(results(1).params.a_min) ', tau-min = ' num2str(results(1).params.tau_min) ', tau-max = ' num2str(results(1).params.tau_max) ', p-spike = ' num2str(results(1).params.p_spike)];
% title(title_string)
% set(gcf, 'Position', [1000 100 1400 900]);

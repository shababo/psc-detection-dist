function plot_scatter_stack(traces, offset_step, varargin)

offset = 0;
stim_start = 1;
time_after_stim = 1; %1000

trial_length = size(traces,2);

for trial = 1:size(traces,1)
    
    offset = offset - 2 * min(traces(trial,:));
    
end

stim_top = 2*max(traces(1,:));
stim_bottom = -offset;

% num_stims = length(find(diff(stims(1,:))))/2;
% change_points = find(diff(stims(1,:)));

% if ~isempty(change_points)
% 
%     trial_length = change_points(end) - change_points(1) + stim_start + time_after_stim;
% 
%     for i = 1:num_stims
% 
%         stim_length = change_points(2*i) - change_points(2*i - 1);
%         this_start = change_points(2*i-1)- change_points(1) + stim_start;
%         rectangle('Position', [this_start stim_bottom stim_length stim_top-stim_bottom],'FaceColor','b','EdgeColor','b')
%         hold on
%     end
% else
%     trial_length = default_length;
% 
% end

offset = 0;


    
for trial = 1:size(traces,1)
    
%     if isempty(change_points)
%         this_trial_start = 1;
%     else
%         this_trial_start = find(stims(trial,:),1,'first') - stim_start;
%     end
    
    trace_to_plot = traces(trial,:);
    sample_inds = find(trace_to_plot);
    scatter(sample_inds/20000,ones(1,length(sample_inds)) - offset +15,(trace_to_plot(sample_inds).^.5)*500,'filled')
    hold on
   
%     if ~isempty(events)
%         scatter((events{trial} - stim_start)/20000,(traces(trial,this_trial_start) - offset - trace_to_plot(1) + offset_step/3)*ones(size(events{trial})),[],colors(ceil(trial/2),:),'filled')
%         hold on
%     end
    
    offset = offset + offset_step;
    
    
end

% axis tight
% axis off

hold off


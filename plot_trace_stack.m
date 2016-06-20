function plot_trace_stack(traces, offset_step, colors, linespec, varargin)


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

if length(varargin) > 1 && ~isempty(varargin{2})
    vert_offset = varargin{2};
else
    vert_offset = 0;
end

if length(varargin) > 2 && ~isempty(varargin{3})
    linewidth = varargin{3};
else
    linewidth = 2;
end
    
for trial = 1:size(traces,1)
    
%     if isempty(change_points)
%         this_trial_start = 1;
%     else
%         this_trial_start = find(stims(trial,:),1,'first') - stim_start;
%     end
    
    trace_to_plot = traces(trial,:);
    %median(trace_to_plot)
    plot((0:trial_length-1),trace_to_plot - offset - trace_to_plot(1) + vert_offset,linespec,'LineWidth',linewidth,'Color',colors(trial,:))
    hold on
   
%     if ~isempty(events)
%         scatter((events{trial} - stim_start)/20000,(traces(trial,this_trial_start) - offset - trace_to_plot(1) + offset_step/3)*ones(size(events{trial})),[],colors(ceil(trial/2),:),'filled')
%         hold on
%     end
    
    offset = offset + offset_step;
    
    
end

if ~isempty(varargin)
    bar_limits = varargin{1};

    if ~isempty(bar_limits)

        bar_corner_time = trial_length/10;
        bar_corner_y = -offset + vert_offset;

        plot([bar_corner_time; bar_corner_time], bar_corner_y + [0; bar_limits(2)], '-k',  bar_corner_time + [0; bar_limits(1)], [bar_corner_y; bar_corner_y], '-k', 'LineWidth', 2)
        text(bar_corner_time - bar_limits(1)/2,bar_corner_y + bar_limits(2)/2, [num2str(bar_limits(2)) ' pA'], 'HorizontalAlignment','right')
        text(bar_corner_time + bar_limits(1)/2,bar_corner_y - bar_limits(2)/2, [num2str(bar_limits(1)*1000) ' ms'], 'HorizontalAlignment','center')
    end
end

axis tight
axis off

hold off


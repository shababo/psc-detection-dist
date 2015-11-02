function plot_scatter_stack(traces, offset_step, bin_edges, varargin)



offset = 0;

if ~isempty(varargin)
    vert_offset = varargin{1};
    if length(varargin) > 1
        multiplier = varargin{2};
    else
        multiplier = 400;
    end
else
    vert_offset = 0;
    multiplier = 400;
end

bin_half_width = .5*(bin_edges(2) - bin_edges(1));

    
for trial = 1:size(traces,1)
    
%     if isempty(change_points)
%         this_trial_start = 1;
%     else
%         this_trial_start = find(stims(trial,:),1,'first') - stim_start;
%     end
    
    traces_to_plot = squeeze(traces(trial,:,:));
    for i = 1:size(traces_to_plot,1);
        trace_to_plot = traces_to_plot(i,:);
        sample_inds = find(trace_to_plot);
        scatter(sample_inds/20000,ones(1,length(sample_inds))*-(bin_edges(i)-bin_half_width) - offset +vert_offset,(trace_to_plot(sample_inds).^.5)*multiplier,'filled')
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


function plot_scatter_stack(traces, offset_step, bin_edges, varargin)



offset = 0;

if ~isempty(varargin) && ~isempty(varargin{1})
    vert_offset = varargin{1};
    if length(varargin) > 1 && ~isempty(varargin{2})
        multiplier = varargin{2};
    else
        multiplier = 400;
    end
else
    vert_offset = 0;
    multiplier = 400;
end

if length(varargin) > 2 && ~isempty(varargin{3})
    colors = varargin{3};
else
    colors = [0 0 1];
end

bin_half_width = .5*(bin_edges(2) - bin_edges(1));

% colors = 1:10;
    
for trial = 1:size(traces,1)
    
%     if isempty(change_points)
%         this_trial_start = 1;
%     else
%         this_trial_start = find(stims(trial,:),1,'first') - stim_start;
%     end
    
    traces_to_plot = squeeze(traces(trial,:,:));
    traces_to_plot = traces_to_plot/max(max(traces_to_plot));
    for i = 1:size(traces_to_plot,1);
        trace_to_plot = traces_to_plot(i,:);
        sample_inds = find(trace_to_plot > .00);
        scatter(sample_inds/20000,ones(1,length(sample_inds))*-(bin_edges(i)-bin_half_width-bin_edges(1)) - offset +vert_offset,ceil(trace_to_plot(sample_inds)*125),colors,'filled')
        hold on
   
%     if ~isempty(events)
%         scatter((events{trial} - stim_start)/20000,(traces(trial,this_trial_start) - offset - trace_to_plot(1) + offset_step/3)*ones(size(events{trial})),[],colors(ceil(trial/2),:),'filled')
%         hold on
%     end
    end
    offset = offset + offset_step;
    
    
end

% axis tight
% axis off

hold off

assignin('base','test',ceil(trace_to_plot(sample_inds)*10))

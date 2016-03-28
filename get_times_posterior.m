function time_posteriors = get_times_posterior(results_file, trace_offset, do_plot, varargin)

load(results_file)
params.traces_filename
try
    
    load(params.traces_filename)
    disp('succesfully loaded')
catch
    disp('loading default')
    load('data/for-paper/direct-stim-w-events-real.mat')

end

% load('/home/shababo/Projects/Mapping/code/psc-detection/data/simulated-data-longer-traces-epsc.mat')
% load('/home/shababo/Desktop/simulated-data-longer-traces-epsc.mat')

if ~isfield(params,'start_ind')
    params.start_ind = 1;
end

if ~isfield(params,'duration')
    params.duration = size(traces,2);
end

if isfield(params,'is_grid') && params.is_grid
    traces = stack_traces(traces);
%     results = stack_results(results);
end

traces = traces(:,params.start_ind:(params.start_ind + params.duration - 1));

if isfield(params,'traces_ind')
    traces = traces(params.traces_ind,:);
end

if length(varargin) > 2 && ~isempty(varargin{3})
    params.traces_ind = varargin{3};
else
    params.traces_ind = 1:size(traces,1);
end
traces = traces(params.traces_ind,:);

if ~isempty(varargin) && ~isempty(varargin{1})
    traces_ind = varargin{1};
    traces = traces(traces_ind,:);
    true_signal = true_signal(traces_ind,:);
end

if length(varargin) >= 2
    max_sample = varargin{2};
end

T = size(traces,2);

if isfield(params,'event_samples')
    event_samples = params.event_samples;
else
    event_samples = 4000;
end

time_posteriors = zeros(size(traces));

burn_in = 1;
% burn_in = 250;


for ii = 1:length(params.traces_ind)
    
    i = params.traces_ind(ii);
    
    if ~exist('max_sample','var')
        map_i = results(i).map_ind;
    else
        [map_val, map_i] = max(results(i).trials.obj(1:max_sample));
    end
    this_posterior = zeros(1,T);
    
    for j = burn_in:length(results(i).trials.times)
        
        these_inds = ceil(results(i).trials.times{j});
        
%         if ~isempty(these_inds)
            this_posterior(these_inds) = this_posterior(these_inds) + 1/length(results(i).trials.times(burn_in:end));
%         end
    end
    
    time_posteriors(ii,:) = this_posterior;
%     for j = 1:length(results(i).trials.times{map_i})
%         
%         
%         ef = genEfilt_ar(results(i).trials.tau{map_i}{j},event_samples);
%         [~, this_curve, ~] = addSpike_ar(results(i).trials.times{map_i}(j),...
%                                             this_curve, 0, ef, ...
%                                             results(i).trials.amp{map_i}(j),...
%                                             results(i).trials.tau{map_i}{j},...
%                                             traces(i,:), ...
%                                             results(i).trials.times{map_i}(j), ...
%                                             2, 1, 1);
%                                         
%     end
%     map_curves(i,:) = this_curve + results(i).trials.base{map_i};
end

if do_plot
    ax1 = axes('Position',[0 0 1 1],'Visible','off');
    ax2 = axes('Position',[.3 .1 .6 .8]);

    descr = {'Parameters:'
        ['a_{min} = ' num2str(params.a_min)];
        ['tau^1_{min} = ' num2str(params.tau1_min)];
        ['tau^1_{max} = ' num2str(params.tau1_max)];
        ['tau^2_{min} = ' num2str(params.tau2_min)];
        ['tau^2_{max} = ' num2str(params.tau2_max)];
        ['p_{spike} = ' num2str(params.p_spike)]
        ['num sweeps = ' num2str(params.num_sweeps)]};


    axes(ax1) % sets ax1 to current axes
    text(.025,0.6,descr)

    axes(ax2)
    plot_trace_stack(traces,trace_offset,bsxfun(@plus,zeros(length(traces),3),[0 0 0]),'-',[.005 10],0)
    hold on
    if exist('true_signal','var')
        times_vec = zeros(size(time_posteriors));
        for i = 1:size(time_posteriors,1)        
            times_vec(i,ceil(true_event_times{i})) = max(max(time_posteriors))+.1;
        end
        plot_scatter_stack(times_vec,trace_offset,[0 0],20,100,[0 0 0])
        hold on
        %plot_trace_stack(true_signal,trace_offset,bsxfun(@plus,zeros(length(traces),3),[0 0 1]),'-',[],80)
        %hold on
    end
    plot_scatter_stack(time_posteriors,trace_offset,[0 0],5,100,[0 0 1])
    hold off

    title(strrep(results_file,'_','-'))

    % if length(varargin) > 1 && varargin{2}
    %     [dir,name,~] = fileparts(results_file);
    %     savefig([dir '/' name '.fig'])
    % end

    % map = params.event_sign*map_curves;
end




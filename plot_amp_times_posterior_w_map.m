function time_amp_posteriors = plot_amp_times_posterior_w_map(results_file, trace_offset, bin_edges, varargin)

load(results_file)

try
    load(params.traces_filename)
catch
    [pathname, filename] = fileparts(params.traces_filename);
    pathname = '/media/shababo/Layover/projects/mapping/code/psc-detection/data/for-paper/';
    load([pathname filename])
end

% load('/home/shababo/Projects/Mapping/code/psc-detection/data/simulated-data-longer-traces-epsc.mat')
% load('/home/shababo/Desktop/simulated-data-longer-traces-epsc.mat')

if ~isfield(params,'start_ind')
    params.start_ind = 1;
end

if ~isfield(params,'duration')
    params.duration = size(traces,2);
end

traces = traces(:,params.start_ind:(params.start_ind + params.duration - 1));

if isfield(params,'traces_ind')
    traces = traces(params.traces_ind,:);
    
end

if ~isempty(varargin) && ~isempty(varargin{1})
    traces_ind = varargin{1};
    traces = traces(traces_ind,:);
    if exist('true_signal','var')
        true_signal = true_signal(traces_ind,:);
    end
    results = results(traces_ind);
end

if length(varargin) >= 2 && ~isempty(varargin{2})
    max_sample = varargin{2};
end


[num_traces, T] = size(traces);
num_amp_bins = length(bin_edges);

if isfield(params,'event_samples')
    event_samples = params.event_samples;
else
    event_samples = 4000;
end

time_amp_posteriors = zeros(size(traces,1),length(bin_edges),T);

for i = 1:size(traces,1)
    
    if ~exist('max_sample','var')
        map_i = results(i).map_ind;
    else
        [map_val, map_i] = max(results(i).trials.obj(1:max_sample));
    end
    this_time_amp_posterior = zeros(length(bin_edges),T);
    
    for j = 1:length(results(i).trials.times)
        
        for k = 1:length(results(i).trials.times{j})
            this_time_ind = ceil(results(i).trials.times{j}(k));
            this_amp_ind = find(results(i).trials.amp{j}(k) < bin_edges,1,'first');
        
            this_time_amp_posterior(this_amp_ind,this_time_ind) =...
                this_time_amp_posterior(this_amp_ind,this_time_ind) + 1/length(results(i).trials.times);
        end
        
    end
    
    time_amp_posteriors(i,:,:) = this_time_amp_posterior;
    this_curve = zeros(1,T);
    for j = 1:length(results(i).trials.times{map_i})
        
%         if results(i).trials.tau{map_i}{j}(2) < 300
        
            ef = genEfilt_ar(results(i).trials.tau{map_i}{j},event_samples);
            [~, this_curve, ~] = addSpike_ar(results(i).trials.times{map_i}(j),...
                                                this_curve, 0, ef, ...
                                                results(i).trials.amp{map_i}(j),...
                                                results(i).trials.tau{map_i}{j},...
                                                traces(i,:), ...
                                                results(i).trials.times{map_i}(j), ...
                                                2, 1, 1);
%         end
                                        
    end
    map_curves(i,:) = this_curve;% + results(i).trials.base{map_i};
end

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
% plot_trace_stack(traces,trace_offset,bsxfun(@plus,zeros(length(traces),3),[1 .4 .4]),'-',[.005 25],0)
% hold on

plot_trace_stack(traces,trace_offset,bsxfun(@plus,zeros(length(traces),3),[0 0 1]),'-',[],100)
hold on
plot_scatter_stack(time_amp_posteriors,trace_offset,bin_edges,0,2000)
hold on
if exist('true_signal','var')
    times_vec = zeros(1,size(time_amp_posteriors,3));
%     times_vec(ceil(true_event_times{1})) = max(max(max(time_amp_posteriors)))+.1;

%     hold on
    plot_trace_stack(true_signal,trace_offset,bsxfun(@plus,zeros(length(traces),3),[0 0 1]),'-',[],0,1)
    hold on
    scatter(true_event_times{1}/20000,-true_amplitudes{1},'*k')
    hold on
end


plot_trace_stack(params.event_sign*map_curves,trace_offset,bsxfun(@plus,zeros(length(traces),3),[1 0 0]),'--',[.01 5],0,1)
hold on



hold off

title(strrep(results_file,'_','-'))

if length(varargin) > 2 && ~isempty(varargin{3})
    [dir,name,~] = fileparts(results_file);
    savefig([dir '/' name '.fig'])
end

% map = params.event_sign*map_curves;
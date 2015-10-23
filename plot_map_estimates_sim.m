function plot_map_estimates_sim(results_file, trace_offset, varargin)

load(results_file)

try
    load(params.traces_filename)
catch
    load('data/sim-temp.mat')
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

map_curves = zeros(size(traces));

for i = 1:size(traces,1)
    
    if ~exist('max_sample','var')
        map_i = results(i).map_ind;
    else
        [map_val, map_i] = max(results(i).trials.obj(1:max_sample));
    end
    this_curve = zeros(1,T);
    
    for j = 1:length(results(i).trials.times{map_i})
        
        ef = genEfilt_ar(results(i).trials.tau{map_i}{j},event_samples);
        [~, this_curve, ~] = addSpike_ar(results(i).trials.times{map_i}(j),...
                                            this_curve, 0, ef, ...
                                            results(i).trials.amp{map_i}(j),...
                                            results(i).trials.tau{map_i}{j},...
                                            traces(i,:), ...
                                            results(i).trials.times{map_i}(j), ...
                                            2, 1, 1);
                                        
    end
    map_curves(i,:) = this_curve + results(i).trials.base{map_i};
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
plot_trace_stack(traces,zeros(size(traces)),[],bsxfun(@plus,zeros(length(traces),3),[1 .4 .4]),[],size(traces,2)-1,trace_offset,'-')
hold on
plot_trace_stack(params.event_sign*map_curves,zeros(size(traces)),[],bsxfun(@plus,zeros(length(traces),3),[0 0 1]),[],size(traces,2)-1,trace_offset,'-')
if exist('true_signal','var')
    hold on
    plot_trace_stack(true_signal,zeros(size(traces)),[],bsxfun(@plus,zeros(length(traces),3),[0 1 0]),[],size(traces,2)-1,trace_offset,'--')
end
hold off

title(strrep(results_file,'_','-'))

if length(varargin) > 1 && varargin{2}
    [dir,name,~] = fileparts(results_file);
    savefig([dir '/' name '.fig'])
end






function plot_map_events(results_file,event_samples,varargin)

load(results_file)

try
    load(params.traces_filename)
catch
    load('data/simulated-epscs-1027.mat')
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

% if isfield(params,'event_samples')
%     event_samples = params.event_samples
% else
%     event_samples = 4000;
% end



% for i = 1:size(traces,1)
    
    if ~exist('max_sample','var')
        map_i = results(1).map_ind;
    else
        [map_val, map_i] = max(results(1).trials.obj(1:max_sample));
    end

    events = zeros(17,event_samples);
    
    
    
    for j = 1:length(results(1).trials.times{map_i})
        
        this_event = zeros(1,event_samples);
        
        ef = genEfilt_ar(results(1).trials.tau{map_i}{j},event_samples);
        [~, this_event, ~] = addSpike_ar([1 1 1 1 1 1],...
                                            this_event, 0, ef, ...
                                            1,...
                                            results(1).trials.tau{map_i}{j},...
                                            traces(1,:), ...
                                            1, ...
                                            2, 1, 1);
        events(j,:) = this_event;
                                        
    end
    
% end

plot(-events','linewidth',2)
axis off

% ax1 = axes('Position',[0 0 1 1],'Visible','off');
% ax2 = axes('Position',[.3 .1 .6 .8]);
% 
% descr = {'Parameters:'
%     ['a_{min} = ' num2str(params.a_min)];
%     ['tau^1_{min} = ' num2str(params.tau1_min)];
%     ['tau^1_{max} = ' num2str(params.tau1_max)];
%     ['tau^2_{min} = ' num2str(params.tau2_min)];
%     ['tau^2_{max} = ' num2str(params.tau2_max)];
%     ['p_{spike} = ' num2str(params.p_spike)]
%     ['num sweeps = ' num2str(params.num_sweeps)]};


% axes(ax1) % sets ax1 to current axes
% text(.025,0.6,descr)

% axes(ax2)
% plot_trace_stack(traces,trace_offset,bsxfun(@plus,zeros(length(traces),3),[1 .4 .4]),'-')
% hold on
% plot_trace_stack(params.event_sign*map_curves,trace_offset,bsxfun(@plus,zeros(length(traces),3),[0 0 1]),'-')
% if exist('true_signal','var')
%     hold on
%     plot_trace_stack(true_signal,trace_offset,bsxfun(@plus,zeros(length(traces),3),[0 1 0]),'--')
% end
% hold off
% 
% title(strrep(results_file,'_','-'))
% 
% if length(varargin) > 1 && varargin{2}
%     [dir,name,~] = fileparts(results_file);
%     savefig([dir '/' name '.fig'])
% end

% map = params.event_sign*map_curves;





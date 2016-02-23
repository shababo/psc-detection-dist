function plot_times_posterior_w_direct_stim(results_file, trace_offset, varargin)

load(results_file)

try
    load(params.traces_filename)
catch e
    disp('loading default traces')
    params.traces_filename = 'data/for-paper/direct-stim-w-events-real.mat';
    load(params.traces_filename)
end

if ~isfield(params,'start_ind')
    params.start_ind = 1;
end

if ~isfield(params,'duration')
    params.duration = size(traces,2);
end

if params.is_grid
    [traces, rebuild_map] = stack_traces(traces);
    results = stack_results(results);
end

traces = traces(:,params.start_ind:(params.start_ind + params.duration - 1));

if isfield(params,'traces_ind')
    traces = traces(params.traces_ind,:);
end

if ~isempty(varargin) && ~isempty(varargin{1})
    traces_ind = varargin{1};
    traces = traces(traces_ind,:);
end

T = size(traces,2);

if isfield(params,'event_samples')
    event_samples = params.event_samples;
else
    event_samples = 4000;
end

map_curves = zeros(size(traces));
mean_curves = zeros(size(traces));
direct_stims = zeros(size(traces));
full_sum = zeros(size(traces));

for i = 1:length(results)
    
    min_i = results(i).map_ind;
    
    
    
%     for k = 1:length(results(i).trials.times)
        this_curve = zeros(1,T);
        k = min_i;

        for j = 1:length(results(i).trials.times{k})

            ef = genEfilt_ar(results(i).trials.tau{k}{j},event_samples);
            [~, this_curve, ~] = addSpike_ar(results(i).trials.times{k}(j),...
                                                this_curve, 0, ef, ...
                                                results(i).trials.amp{k}(j),...
                                                results(i).trials.tau{k}{j},...
                                                traces(i,:), ...
                                                results(i).trials.times{k}(j), ...
                                                2, 1, 1);

        end

        if k == min_i
            map_curves(i,:) = params.event_sign*(this_curve + results(i).trials.base{min_i});
        end
        
        mean_curves(i,:) = mean_curves(i,:) - this_curve/length(results(i).trials.times);
%     end

    
    nBins_stim = length(params.stim_in);
    nBins = size(traces,2);
    if nBins_stim < nBins
        stim_in = [params.stim_in zeros(1,nBins - length(params.stim_in))];
    else
        stim_in = params.stim_in(1:nBins);
    end
    t = 0:nBins-1;
    results(i).trials.stim_tau_fall{min_i}
    results(i).trials.stim_tau_rise{min_i}
    stim_decay = exp(-t/results(i).trials.stim_tau_fall{min_i});
    stim_rise = -exp(-t/results(i).trials.stim_tau_rise{min_i});
    stim_kernel = (stim_decay + stim_rise)/sum(stim_decay + stim_rise);
    stim_template = conv(stim_in,stim_kernel);
    stim_template = stim_template(1:nBins)/max(stim_template(1:nBins));
    
    direct_stims(i,:) = -stim_template*results(i).trials.stim_amp{min_i} + results(i).trials.base{min_i};
    
    full_sum(i,:) = map_curves(i,:) + direct_stims(i,:) - results(i).trials.base{min_i};

    
%     results(i).trials.stim_amp{min_i}
end

time_posteriors = get_times_posterior(results_file, 0, 0);

early_cut = 20;

for i = 1:size(time_posteriors,1)
    time_posteriors(i,:) = time_posteriors(i,:)/max(time_posteriors(i,early_cut+1:end));
    time_posteriors(i,1:early_cut) = 0;
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

colors = lines(85);


axes(ax2)
plot_trace_stack(traces,trace_offset,bsxfun(@plus,zeros(length(traces),3),[0 0 0]),'-',[.005 25])
hold on
% plot_trace_stack(full_sum,trace_offset,bsxfun(@plus,zeros(length(direct_stims),3),[.75 0 0]),'--',[],[],2)
hold on
the_color = 42;
plot_trace_stack(direct_stims,trace_offset,bsxfun(@plus,zeros(length(direct_stims),3),colors(the_color,:)),'-',[],25,2)
hold on
the_color = 43;
plot_trace_stack(map_curves,trace_offset,bsxfun(@plus,zeros(length(direct_stims),3),colors(the_color,:)),'-',[],55,2)
hold on
plot_trace_stack(time_posteriors*30,trace_offset,bsxfun(@plus,zeros(length(traces),3),[.1 .1 .75]),'-',[],80,1)
hold off

title(strrep(results_file,'_','-'))






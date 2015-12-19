%%
clc; clear all
rng(8423)
load('data/good-example-trace-1000.mat')
 load('/media/shababo/Layover/projects/mapping/code/psc-detection/data/for-paper/epscs/L5andL23PC_data_spontaneous_EPSCs.mat')
phi = results(1).trials.phi{results(1).map_ind}
sigmasq = results(1).trials.noise{results(1).map_ind};
simulate_data; sim_ar0_data = traces;

true_signal_tmp(1,:) = zeros(size(true_signal));
true_event_times_tmp = cell(3,1);
true_amplitudes_tmp = cell(3,1);
true_taus_tmp = cell(3,1);

true_signal_tmp(2,:) = true_signal;
true_event_times_tmp{2} = true_event_times{1};
true_amplitudes_tmp{2} = amplitudes{1}; true_taus_tmp{2} = taus{1};

load('data/good-example-trace-0009.mat')



phi = results(1).trials.phi{results(1).map_ind};
sigmasq = results(1).trials.noise{results(1).map_ind};
simulate_data; sim_ar2_data = traces;
% load('data/for-paper/good-example-trace.mat')
i = 12; start_ind = max(1,ceil(rand*length(Master(4).spontSweep{i}))-3050);
traces = [Master(4).spontSweep{i}(start_ind:(start_ind+2999))'; sim_ar0_data; sim_ar2_data];
figure;

plot_trace_stack(traces,30,zeros(3,3),'-',[.010 10])

true_signal_tmp(3,:) = true_signal;
true_event_times_tmp{3} = true_event_times{1};
true_amplitudes_tmp{3} = amplitudes{1}; true_taus_tmp{3} = taus{1};

true_signal = true_signal_tmp;
true_event_times = true_event_times_tmp;
true_amplitudes = true_amplitudes_tmp;
true_taus = true_taus_tmp;
%%

% load('data/simulated-epscs-1027.mat')
template = zeros(1,600);
t = 0:599;
num_events = 0;

for j = 1:length(true_amplitudes)
    j
    for i = 1:length(true_amplitudes{j})
        i
        tau_decay = true_taus{j}{i}(2); decay = exp(-t/tau_decay);
        tau_rise = true_taus{j}{i}(1); rise = -exp(-t/tau_rise);
        % plot(decay); hold on; plot(rise)
        template = template + -(decay + rise)/max(decay+rise);
        num_events = num_events + 1;
    end
end

template = template/num_events;
%%
save('data/harder-for-cosyne.mat','traces','true_signal','true_event_times','true_amplitudes','true_taus','template')



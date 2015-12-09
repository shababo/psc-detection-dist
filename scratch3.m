alex_l23_pyr = get_sweeps('/home/shababo/projects/mapping/code/psc-detection/data/for-paper/epscs/L5andL23PC_data_spontaneous/20141016/20141016ExpA.mat',1,1,0);
som_0612_trace = get_sweeps('/home/shababo/projects/mapping/code/psc-detection/data/for-paper/epscs/som/0612_som01.mat',1,3,0);
alan_emx = get_sweeps('/home/shababo/projects/mapping/code/psc-detection/data/for-paper/ipscs/EMX/20150706_AM2F_TI1_cell6-7.mat',1,5,0);

figure;
subplot(311)
plot(alex_l23_pyr)
subplot(312)
plot(som_0612_trace);
subplot(313)
plot(alan_emx)

%%
load('data/for-paper/real_noise_traces_longer.mat')
experiment = '0000';
load(['longer-noise-examples-epscs-' experiment '.mat']);


phi_ar2_noise_01 = results(1).trials.phi{results(1).map_ind};
phi_ar2_noise_02 = results(2).trials.phi{results(2).map_ind};
sigmasq_ar2_noise_01 = results(1).trials.noise{results(1).map_ind};
sigmasq_ar2_noise_02 = results(2).trials.noise{results(2).map_ind};

load(['longer-noise-examples-ipscs-' experiment '.mat']);

phi_ar2_noise_03 = results(1).trials.phi{results(1).map_ind};
sigmasq_ar2_noise_03 = results(1).trials.noise{results(1).map_ind};
sim_data = zeros(3,2000);
phi = phi_ar2_noise_01; sigmasq = sigmasq_ar2_noise_01; simulate_data; sim_data(1,:) = traces;
phi = phi_ar2_noise_02; sigmasq = sigmasq_ar2_noise_02; simulate_data; sim_data(2,:) = traces;
phi = phi_ar2_noise_03; sigmasq = sigmasq_ar2_noise_03; simulate_data; sim_data(3,:) = traces;
%%
experiment = '1000';
load(['longer-noise-examples-epscs-' experiment '.mat']);


phi_ar0_noise_01 = results(1).trials.phi{results(1).map_ind};
phi_ar0_noise_02 = results(2).trials.phi{results(2).map_ind};
sigmasq_ar0_noise_01 = results(1).trials.noise{results(1).map_ind};
sigmasq_ar0_noise_02 = results(2).trials.noise{results(2).map_ind};

load(['longer-noise-examples-ipscs-' experiment '.mat']);

phi_ar0_noise_03 = results(1).trials.phi{results(1).map_ind};
sigmasq_ar0_noise_03 = results(1).trials.noise{results(1).map_ind};
sim_ar0_data = zeros(3,2000);
phi = phi_ar0_noise_01; sigmasq = sigmasq_ar0_noise_01; simulate_data; sim_ar0_data(1,:) = traces;
phi = phi_ar0_noise_02; sigmasq = sigmasq_ar0_noise_02; simulate_data; sim_ar0_data(2,:) = traces;
phi = phi_ar0_noise_03; sigmasq = sigmasq_ar0_noise_03; simulate_data; sim_ar0_data(3,:) = traces;

%%
figure;
subplot(131)
plot_trace_stack(real_noise_traces_longer(:,6001:8000),20,zeros(3,3),'-',[.010 10])
title('Voltage Clamp Recordings')

subplot(132)
plot_trace_stack(sim_ar0_data,20,zeros(3,3),'-',[.010 10])
title('Simulated Noise From Fits - AR(0)')

subplot(133)
plot_trace_stack(sim_data,20,zeros(3,3),'-',[.010 10])
title('Simulated Noise From Fits - AR(2)')

%%
rng(18480)
load('data/good-example-trace-1000.mat')


phi = results(1).trials.phi{results(1).map_ind};
sigmasq = results(1).trials.noise{results(1).map_ind};
simulate_data; sim_ar0_data = traces;

load('data/good-example-trace-0009.mat')


phi = results(1).trials.phi{results(1).map_ind};
sigmasq = results(1).trials.noise{results(1).map_ind};
simulate_data; sim_ar2_data = traces;
% load('data/for-paper/good-example-trace.mat')
 load('data/for-paper/epscs/L5andL23PC_data_spontaneous_EPSCs.mat')

traces_to_plot = [Master(4).spontSweep{13}(1028:3027)'; sim_ar0_data; sim_ar2_data];

figure;
plot_trace_stack(traces_to_plot,20,zeros(3,3),'-',[.010 10])

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


template = zeros(1,600);
t = 0:599;
num_events = 0;
for j = 2:length(true_amplitudes)
    for i = 1:length(true_amplitudes{j})

        tau_decay = true_taus{j}{i}(2); decay = exp(-t/tau_decay);
        tau_rise = true_taus{j}{i}(1); rise = -exp(-t/tau_rise);
        % plot(decay); hold on; plot(rise)
        template = template + -(decay + rise)/max(decay+rise);
        num_events = num_events + 1;
    end
end

template = template/num_events;

% save('data/for-paper/real-vs-ar0-vs-ar2-sim-cosyne-abs.mat','traces','true_signal','true_event_times','true_amplitudes','true_taus','template')




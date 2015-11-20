%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate data from model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rng(12341)

% Simulate K voltage-clamp observed neurons
K = 1;

% Simulation: Given sampling, what indicator timecourse. Also depends on spike rate
% Sampling rate
% Noise level
% Time constants
% Firing rate
% Poisson or periodic

T = 4000; %bins - start not too long
binSize = 1/20000; %
tau_r_bounds = [5 20];
tau_f_bounds = [20 150];
firing_rate = 0; %spike/sec 
% c_noise = 2.5; %calcium signal std
baseline = 0;
A = 1; %magnitude scale of spike response

Y = zeros(K,T);
C = zeros(K,T);
Y_AR = zeros(K,T);
Spk = cell(1,K);
taus = cell(1,K);
amplitudes = cell(1,K);

periodic = 0; %if zero, uses poisson spiketrain.

tau = [mean(tau_r_bounds)/binSize mean(tau_f_bounds)/binSize]; %time constants in bin-units

% compute exponential filter(s)
ef=genEfilt(tau,T);

n_spike = firing_rate*(T*binSize);
p_spike = n_spike/T;

times = cumsum(binSize*1e-3*ones(1,T),2); % in sec

a_min = 2.5;
a_max = 8;
nc = 1; %trials

user_defined = 1;

for ki = 1:K
    
%     if ~exist('true_signal','var')
        ssi = [];
        % startingSpikeTimes = [10 20];
        startingSpikeTimes = [];
        if periodic
            s_int = T/n_spike;
            startingSpikeTimes = [startingSpikeTimes s_int:s_int:(T-10)];
        elseif user_defined
            startingSpikeTimes = 200:400:T;
            amps = (1:length(startingSpikeTimes))*.5;
        else
            num_spikes = poissrnd(n_spike);
            startingSpikeTimes = T*rand(1,num_spikes);
    %         startingSpikeTimes = times(rand(1,T)<p_spike)/(binSize*1e-3);
        end
        ci = baseline*ones(nc,T); %initial calcium is set to baseline 

        offsets = zeros(nc,1);
        rescalings = ones(nc,1);

        st_std = .1; %across trials
        ati = cell(nc,1); % array of lists of individual trial spike times
        ati_ = cell(nc,1); 
        sti = cell(nc,1); 
        sti_ = cell(nc,1); 

        N = 0;
        Dt = 1;

        trace_amps = [];
        trace_taus = {};

        for i = 1:length(startingSpikeTimes)        
            tmpi = startingSpikeTimes(i); 
            ssi_ = [ssi tmpi];
            cti_ = ci;

            ati_ = ati;
            logC_ = 0;
            for ti = 1:nc
                a_init = amps(i);
                tmpi_ = tmpi+(st_std*randn);
%                 tau(1) = diff(tau_r_bounds)*rand() + tau_r_bounds(1);
%                 tau(2) = diff(tau_f_bounds)*rand() + tau_f_bounds(1);
                tau(1) = 15;
                tau(2) = 100;
                ef=genEfilt(tau,T);
                [si_, ci_, logC_] = addSpike(sti{ti},ci(ti,:),logC_,ef,a_init,tau,ci(ti,:),tmpi_, N+1, Dt, A); %adds all trials' spikes at same time
                sti_{ti} = si_;
                cti_(ti,:) = ci_;
                ati_{ti} = [ati_{ti} a_init];
            end
            ati = ati_;
            ssi = ssi_;
            sti = sti_;
            ci = cti_;
            logC = logC_;
            N = length(ssi); %number of spikes in spiketrain
            trace_amps = [trace_amps a_init];
            trace_taus{i} = tau;
        end

        % add direct stim - SET TO ZERO RIGHT NOW
        stim_tau_rise = .0015*20000; % values for chr2 from lin et al 2009 (biophysics)
        stim_tau_fall = .013*20000;
        stim_amp = 0;
        stim_start = 500;
        stim_duration = .05*20000;
        stim_in = [zeros(1,stim_start) ones(1,stim_duration) zeros(1,T-stim_start-stim_duration)];
        t = 0:T-1;
        stim_decay = exp(-t/stim_tau_fall);
        stim_rise = -exp(-t/stim_tau_rise);
        stim_kernel = (stim_decay + stim_rise)/sum(stim_decay + stim_rise);
        stim_response = conv(stim_in,stim_kernel);
        stim_response = stim_amp*stim_response(1:T)/max(stim_response(1:T));
        d = ci + stim_response;
%     else
%         disp('USING PREVIOUSLY MADE SIGNAL')
%     end

%     sigmasq = 2.0;
    c_noise = sqrt(sigmasq);
    
%     phi = [1, 1.0, -.42]; %this determines what the AR noise looks like.
    p = length(phi) - 1;
    U = c_noise*randn(nc,T);
    er = zeros(T,1);

    for t = (1+p):(T+p)
        er(t) = phi*[U(t-p); er(t-1:-1:(t-p))];
    end

    er = er((1+p):(T+p))';
    y = d + U;
    y_ar = d + er;

    
    Y(ki,:) = y;
    C(ki,:) = ci;
    Y_AR(ki,:) = y_ar;
    Spk{ki} = startingSpikeTimes;
    amplitudes{ki} = trace_amps;
    taus{ki} = trace_taus;
end

figure;
ax1 = subplot(311);
plot_trace_stack((-C' - 20*repmat(0:(K-1),T,1))',0,zeros(K,3),'-',[])
title('True Current')

ax2 = subplot(312);
% plot((0:T-1)/20000,-Y' - 20*repmat(0:(K-1),T,1))
plot_trace_stack((-er' - 20*repmat(0:(K-1),T,1))',0,zeros(K,3),'-',[])
title('AR(2) Noise Process')

ax3 = subplot(313);
% plot((0:T-1)/20000,-Y_AR' - 20*repmat(0:(K-1),T,1))
plot_trace_stack((-Y_AR' - 20*repmat(0:(K-1),T,1))',25,zeros(K,3),'-',[.01 10])
title(['Observation'])
xlimits = get(gca,'xlim');
ylimits = get(gca,'ylim');
set(ax2,'xlim',xlimits)
set(ax2,'ylim',ylimits)
set(ax1,'xlim',xlimits)
set(ax1,'ylim',ylimits)




% xlim([1 2000])
% ylim([-20 20])

%ax4 = subplot(313);
%plot(map_est)
% xlim([1 2000])
% ylim([-20 20])

% linkaxes([ax1 ax2 ax3])

traces = -Y_AR;
true_signal = -C;
true_event_times = Spk;
true_amplitudes = amplitudes; true_taus = taus;

sorted_times = sort(true_event_times{1});
legend_names = {};
%%
% figure;
% % subplot(211)
% for i = 1:length(true_amplitudes{1})
%     t = 0:1:.015*20000;
%     tau_decay = taus{1}{i}(2); decay = exp(-t/tau_decay);
%     tau_rise = taus{1}{i}(1); rise = -exp(-t/tau_rise);
%     % plot(decay); hold on; plot(rise)
%     plot(-(decay + rise)/max(decay+rise)*true_amplitudes{1}(i),'linewidth',2)
%     hold on; 
%     legend_names{i} = ['event ' num2str(find(sorted_times == true_event_times{1}(i)))];
% end
% hold on
% legend(legend_names)
% 
% bar_limits = [.0025*20000 2.5];
% bar_corner_time = t(end)/6;
% bar_corner_y = 0;
% offset = 15;
% 
% plot([bar_corner_time; bar_corner_time], -offset + [0; bar_limits(2)], '-k',  bar_corner_time + [0; bar_limits(1)], [-offset; -offset], '-k', 'LineWidth', 2)
% text(bar_corner_time - bar_limits(1)/2,bar_corner_y + bar_limits(2)/2 - offset, [num2str(bar_limits(2)) ' pA'], 'HorizontalAlignment','right')
% text(bar_corner_time + bar_limits(1)/2,bar_corner_y - bar_limits(2)/2 - offset, [num2str(bar_limits(1)/20000*1000) ' ms'], 'HorizontalAlignment','center')
% 
% axis off
%%
%     % subplot(212)
%     % t = 0:1:.02*20000;
%     tau_decay = taus{1}{2}(2); decay = exp(-t/tau_decay);
%     tau_rise = taus{1}{2}(1); rise = -exp(-t/tau_rise);
%     % plot(decay); hold on; plot(rise)
%     hold on; plot((decay + rise)/max(decay+rise)*true_amplitudes{1}(2))


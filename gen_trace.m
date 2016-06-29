function [trace, true_signal] = gen_trace(data_params,bg_params,evoked_params)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate data from model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rng(12341)

% Simulate K voltage-clamp observed neurons

K = 1;

% Simulation: Given sampling, what indicator timecoparamsurse. Also depends on spike rate
% Sampling rate
% Noise level
% Time constants
% Firing rate
% Poisson or periodic

% data_params.data_params.T = 2000; %bins - start not too long
% data_params.data_params.dt = 1/20000; %
% data_params.data_params.baseline = 0;
% data_params.sigmasq = 3.5;
% data_params.data_params.phi = [1, .80, -.12]; %this determines what the AR noise looks like.
c_noise = sqrt(data_params.sigmasq);
% 
% 
% bg_params.bg_params.tau_r_bounds = [5 20];
% bg_params.bg_params.tau_f_bounds = [20 150];
% bg_params.bg_params.a_min = .5;
% bg_params.bg_params.a_max = 6;
% bg_params.bg_params.firing_rate = 0; %spike/sec 
% 
% evoked_params.evoked_params.stim_tau_rise = .0015*20000; % values for chr2 from lin et al 2009 (biophysics)
% evoked_params.evoked_params.stim_tau_fall = .013*20000;
% evoked_params.evoked_params.stim_amp = 0;
% evoked_params.evoked_params.stim_start = .05*20000;
% evoked_params.evoked_params.stim_duration = .05*20000;
% evoked_params.times
% evoked_params.a
% evoked_params.tau_r
% evoked_params.tau_f

A = 1; %magnitude scale of spike response

Y = zeros(1,data_params.T);
C = zeros(1,data_params.T);
Y_AR = zeros(1,data_params.T);
Spk = cell(1,K);
taus = cell(1,K);
amplitudes = cell(1,K);


tau = [mean(bg_params.tau_r_bounds)/data_params.dt mean(bg_params.tau_f_bounds)/data_params.dt]; %time constants in bin-units

% compute exponential filter(s)
ef=genEfilt(tau,data_params.T);

n_spike = bg_params.firing_rate*(data_params.T*data_params.dt);
p_spike = n_spike/data_params.T;

times = cumsum(data_params.dt*1e-3*ones(1,data_params.T),2); % in sec



nc = 1; %trials


    
%     if ~exist('true_signal','var')
        ssi = [];
        % startingSpikeTimes = [10 20];
        startingSpikeTimes = [];

    num_spikes = poissrnd(n_spike);
    startingSpikeTimes = data_params.T*rand(1,num_spikes);
%         startingSpikeTimes = times(rand(1,data_params.T)<p_spike)/(data_params.dt*1e-3);

    ci = data_params.baseline*ones(nc,data_params.T); %initial calcium is set to data_params.baseline 


    st_std = 0; %across trials
    ati = cell(nc,1); % array of lists of individual trial spike times
    ati_ = cell(nc,1); 
    sti = cell(nc,1); 
    sti_ = cell(nc,1); 

    N = 0;
%         Dt = 1;

    trace_amps = [];
    trace_taus = {};

    % generate bg events
    for i = 1:length(startingSpikeTimes)        
        tmpi = startingSpikeTimes(i); 
        ssi_ = [ssi tmpi];
        cti_ = ci;

        ati_ = ati;
        logC_ = 0;
        for ti = 1:nc
            a_init = bg_params.a_min + (bg_params.a_max-bg_params.a_min)*rand;
            tmpi_ = tmpi+(st_std*randn);
            tau(1) = diff(bg_params.tau_r_bounds)*rand() + bg_params.tau_r_bounds(1);
            tau(2) = diff(bg_params.tau_f_bounds)*rand() + bg_params.tau_f_bounds(1);
            ef=genEfilt(tau,data_params.T);
            [si_, ci_, logC_] = addSpike(sti{ti},ci(ti,:),logC_,ef,a_init,tau,ci(ti,:),tmpi_, N+1, 1, A); %adds all trials' spikes at same time
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
    
    % evoked spikes
    ssi = [];
    ati = cell(nc,1); % array of lists of individual trial spike times
    ati_ = cell(nc,1); 
    sti = cell(nc,1); 
    sti_ = cell(nc,1); 
    for i = 1:length(evoked_params.times)        
        tmpi = evoked_params.times(i); 
        ssi_ = [ssi tmpi];
        cti_ = ci;

        ati_ = ati;
        logC_ = 0;
        for ti = 1:nc
            a_init = evoked_params.a(i);
%             tmpi_ = tmpi+(st_std*randn);
            tau(1) = evoked_params.tau_r(i);
            tau(2) = evoked_params.tau_f(i);
            ef=genEfilt(tau,data_params.T);
            [si_, ci_, logC_] = addSpike(sti{ti},ci(ti,:),logC_,ef,a_init,tau,ci(ti,:),tmpi_, N+1, 1, A); %adds all trials' spikes at same time
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
    stim_in = [zeros(1,evoked_params.stim_start) ones(1,evoked_params.stim_duration) zeros(1,data_params.T-evoked_params.stim_start-evoked_params.stim_duration)];
    t = 0:data_params.T-1;
    stim_decay = exp(-t/evoked_params.stim_tau_fall);
    stim_rise = -exp(-t/evoked_params.stim_tau_rise);
    stim_kernel = (stim_decay + stim_rise)/sum(stim_decay + stim_rise);
    stim_response = conv(stim_in,stim_kernel);
    stim_response = evoked_params.stim_amp*stim_response(1:data_params.T)/max(stim_response(1:data_params.T));
    d = ci + stim_response;
%     else
%         disp('USING PREVIOUSLY MADE SIGNAL')
%     end


p = length(data_params.phi) - 1;
U = c_noise*randn(nc,data_params.T);
er = zeros(data_params.T,1);

for t = (1+p):(data_params.T+p)
    er(t) = data_params.phi*[U(t-p); er(t-1:-1:(t-p))];
end

er = er((1+p):(data_params.T+p))';
y = d + U;
y_ar = d + er;


Y = y;
C = ci;
Y_AR = y_ar;
Spk = startingSpikeTimes;



% figure;
% ax1 = subplot(311);
% plot_trace_stack((-C' - 20*repmat(0:(K-1),data_params.T,1))',0,zeros(K,3),'-',[])
% title('True Current')
% 
% ax2 = subplot(312);
% % plot((0:data_params.T-1)/20000,-Y' - 20*repmat(0:(K-1),data_params.T,1))
% plot_trace_stack((-er' - 20*repmat(0:(K-1),data_params.T,1))',0,zeros(K,3),'-',[])
% title('AR(2) Noise Process')
% 
% ax3 = subplot(313);
% % plot((0:data_params.T-1)/20000,-Y_AR' - 20*repmat(0:(K-1),data_params.T,1))
% plot_trace_stack((-Y_AR' - 20*repmat(0:(K-1),data_params.T,1))',25,zeros(K,3),'-',[.01 10])
% title(['Observation'])
% xlimits = get(gca,'xlim');
% ylimits = get(gca,'ylim');
% set(ax2,'xlim',xlimits)
% set(ax2,'ylim',ylimits)
% set(ax1,'xlim',xlimits)
% set(ax1,'ylim',ylimits)




% xlim([1 2000])
% ylim([-20 20])

%ax4 = subplot(313);
%plot(map_est)
% xlim([1 2000])
% ylim([-20 20])

% linkaxes([ax1 ax2 ax3])

trace = -Y_AR;
true_signal.trace = -C;
true_signal.event_times = Spk;
true_signal.amplitudes = amplitudes; true_signal.taus = taus;
% 
% sorted_times = sort(true_event_times{1});
% legend_names = {};
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


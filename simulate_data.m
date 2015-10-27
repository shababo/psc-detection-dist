%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate data from model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(1234)

% Simulate K calcium imaging observed neurons
K = 1;

% Simulation: Given sampling, what indicator timecourse. Also depends on spike rate
% Sampling rate
% Noise level
% Time constants
% Firing rate
% Poisson or periodic

T = 20000; %bins - start not too long
binSize = 1/20000; %ms
tau_r_bounds = [1 10];
tau_f_bounds = [10 100];
firing_rate = 10; %spike/sec 
c_noise = 2.0; %calcium signal std
baseline = 0;
A = 1; %magnitude scale of spike response

Y = zeros(K,T);
C = zeros(K,T);
Y_AR = zeros(K,T);
Spk = cell(1,K);
taus = cell(1,K);
amplitudes = cell(1,K);

periodic = 0; %if zero, uses poisson spiketrain.

tau = [mean(tau_r_bounds)/binSize mean(tau_f_bounds)/binSize]; %time constants in bin-uni20ts

% compute exponential filter(s)
ef=genEfilt(tau,T);

n_spike = firing_rate*(T*binSize);
p_spike = n_spike/T;

times = cumsum(binSize*1e-3*ones(1,T),2); % in sec

a_min = 2;
a_max = 15;
nc = 1; %trials

for ki = 1:K
    ssi = [];
    % startingSpikeTimes = [10 20];
    startingSpikeTimes = [];
    if periodic
        s_int = T/n_spike;
        startingSpikeTimes = [startingSpikeTimes s_int:s_int:(T-10)];
    else
        startingSpikeTimes = times(rand(1,T)<p_spike)/(binSize*1e-3);
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
            a_init = a_min + (a_max-a_min)*rand;
            tmpi_ = tmpi+(st_std*randn);
            tau(1) = diff(tau_r_bounds)*rand() + tau_r_bounds(1);
            tau(2) = diff(tau_f_bounds)*rand() + tau_f_bounds(1);
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
    
    % add direct stim
    stim_tau_rise = 5;
    stim_tau_fall = 400;
    stim_amp = 50;
    stim_start = 5000;
    stim_duration = 1000;
    stim_in = [zeros(1,stim_start) ones(1,stim_duration) zeros(1,T-stim_start-stim_duration)];
    t = 0:T-1;
    stim_decay = exp(-t/stim_tau_fall);
    stim_rise = -exp(-t/stim_tau_rise);
    stim_kernel = (stim_decay + stim_rise)/sum(stim_decay + stim_rise);
    stim_response = stim_amp*conv(stim_in,stim_kernel);
    stim_response = stim_response(1:T);
    d = ci + stim_response;


    
    c_noise = 2.0;
    p = 2;
    phi = [1, .3, .35]; %this determines what the AR noise looks like.
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
plot(-C' - 10*repmat(0:(K-1),T,1))

ax2 = subplot(312);
% plot(-Y' - 20*repmat(0:(K-1),T,1))

% ax3 = subplot(413);
plot(-Y_AR' - 20*repmat(0:(K-1),T,1))
% xlim([1 2000])
% ylim([-20 20])

ax4 = subplot(313);
plot(map_est)
% xlim([1 2000])
% ylim([-20 20])

% linkaxes([ax2 ax3])

traces = -Y_AR;


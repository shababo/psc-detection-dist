function [trace, true_signal] = sample_traces(K, params)
% SAMPLE_TRACE sample a random voltage-clamp trace
%   [trace, true_signal] = SAMPLE_TRACE(N, params) returns K voltage-clamp
%   traces drawn from the distribution defined by params.
%   params is a struct with these fields:
%       params.T                length of trace in samples
%       params.dt               sample duration in seconds
%       params.phi              a 1 x p vector of temporal dependencies for the
%                               AR(p) noise model
%       params.sigma_sq         variance of noise (pA^2 - or trace unit^2)
%       params.baseline_bounds  a 2 x 1 vector for the bounds of the
%                               uniform prior on the basline (holding
%                               current, in pA - or trace unit)
%       params.tau_r_bounds     a 2 x 1 vector for the bounds of the
%                               uniform prior on rise time constants (in
%                               secs)
%       params.tau_f_bounds     a 2 x 1 vector for the bounds of the
%                               uniform prior on fall time constants (in
%                               secs)
%       params.a_bounds         a 2 x 1 vector for the bounds of the
%                               uniform prior on amplitudes (pA - or trace unit)
%       params.rate              scalar for the Poisson event rate
%                               (events/sec)
%       params.event_direction  1 for IPSCs and -1 for EPSCs, corresponding
%                               to outward and inward currents respectively 


% matrix of noisy data
Y = zeros(K,params.T);

event_times = cell(K,1);
taus = cell(K,1);
amplitudes = cell(K,1);

% compute the expected number of events
e_num_events = params.rate*(params.T*params.dt);
% % draw number of events for each trace
% num_events = poissrnd(e_num_events,K);


% draw baselines
baseline = unifrnd(params.baseline_bounds(1),params.baseline_bounds(2),K,1);
% initialize noiseless traces to baselines
C = bsxfun(@plus, zeros(K,params.T), baseline);  

% number of timesteps for noise filtering
p = length(params.phi);

% generate each trace
for k = 1:K
    
    inter_event_times = -log(rand(2*e_num_events,1))/params.rate/params.dt;
    these_event_times = cumsum(inter_event_times);
    these_event_times = these_event_times(these_event_times < params.T);
    
    % we need to build the times vector as we go because of the way the add
    % event function is built - which is designed more with the Gibbs
    % sampler in mind
    N = 0;
    these_times_tmp = [];
    
    % init for amps and taus
    these_amps = zeros(length(these_event_times),1);
    these_taus = zeros(length(these_event_times),2);
    
    % generate each event
    for i = 1:length(these_event_times) 
        
        % again - this is bookeeping related to how add_event(...) works
        this_time = these_event_times(i); 
              
        % draw event amplitude
        this_amp = params.a_bounds(1) + (params.a_bounds(2)-params.a_bounds(1))*rand;
        
        % draw this event taus
        this_tau(1) = diff(params.tau_r_bounds/params.dt)*rand() + params.tau_r_bounds(1)/params.dt;
        this_tau(2) = diff(params.tau_f_bounds/params.dt)*rand() + params.tau_f_bounds(1)/params.dt;
        
        % gen exp filter
        ef=genEfilt(this_tau,params.T);
        
        % add event to trace
        [these_times_tmp, C(k,:)] = add_event(these_times_tmp, C(k,:), 0, ef, this_amp, this_tau, C(k,:), this_time, N+1);
        
        % update holder vecs
        N = N + 1;
        these_amps(i) = this_amp;
        these_taus(i,:) = this_tau;

    end
    
    % generate white gaussian noise
    white_noise = sqrt(params.sigma_sq)*randn(1,params.T);
    % filter to create ar(p) noise
    ar_noise = zeros(1,params.T+p);
    for t = (1+p):(params.T+p)
        ar_noise(t) = [1 params.phi]*[white_noise(t-p); ar_noise(t-1:-1:(t-p))'];
    end
    ar_noise = ar_noise((1+p):(params.T+p));
    
    Y(k,:) = C(k,:) + ar_noise;
  
    event_times{k} = these_event_times*params.dt;
    taus{k} = these_taus*params.dt;
    amplitudes{k} = these_amps;
    
end

% setup output structs
trace = params.event_direction * Y;
true_signal.trace = params.event_direction * C;
true_signal.event_times = event_times;
true_signal.amplitudes = amplitudes; 
true_signal.taus = taus;



function [posterior, mcmc]  = sample_params(trace, params, times_init, amps_init, taus_init)
% SAMPLE_PARAMS Uses the Gibbs sampler described in Merel et al, J. Neuro.
% Meths. 2016 to sample from the posterior of the latent variables given
% some noisy physiological trace (e.g. a voltage-clamp or trace imaging
% trace)
%   [posterior, mcmc] = SAMPLE_PARAMS(trace, params, times_init, amps_init,
%   taus_init) samples the posterior distribution of the parameters given
%   some phsyiological trace. The inputs to the function are
%       trace       a 1 x T vector of the noisy data
%       params      a struct containging all of the necessary values for
%                   running the algorithm. See the code for GET_PARAMS for
%                   details.
%       times_init  an N x 1 vector of the initialization times
%       amps_init   an N x 1 vector of the initialization amplitudes
%       taus_init   an N x 1 vector of the initialization time constants
%   where K is the number of initial events.
%
%   The outputs from this function are posterior and mcmc. The struct
%   posterior has these fields:
%       amp         a 1 x L* vector of amplitude samples
%       base        a 1 x L vector of baseline samples
%       tau1        a 1 x L* vector of rise tau samples
%       tau2        a 1 x L* vector of fall tau samples
%       num_events  a 1 x L vector of the number of event samples
%       phi         an L x 3 matrix of phi samples
%       noise       a 1 x L vector of samples for the noise variance
%       obj         a 1 x L vector of the objective function value for each sweep
%       times       a 1 x L* vector of the event time samples
%   where L is the number of Gibbs sweeps (including burn-in) and L* is the
%   the sum number of infered events for each sweep (i.e.
%   sum(posterior.num_events). This structure is more memory efficient than
%   other options. Using this structure is not too difficult. The features
%   postiorior.amp(i), postiorior.tau1(i), postiorior.tau2(i), and 
%   postiorior.times(i) all come from the same event. And the events from
%   i-th sweep can be indexed like this (ignoring boundary cases):
%       posterior.amp(sum(posterior.num_events(1:i-1))+1:...
%                  sum(posterior.num_events(1:i-1))+posterior.num_events(i))
%   The output mcmc is a stuct which contains some information on the
%   accept/reject statistics for the sampler. It has fields:
%       addMoves    a 1 x 2 vector where the first entry is the number of
%                   accepts and the second is the total number of proposals
%                   for the add proposals
%       timeMoves    a 1 x 2 vector where the first entry is the number of
%                   accepts and the second is the total number of proposals
%                   for the time proposals
%       dropMoves   a 1 x 2 vector where the first entry is the number of
%                   accepts and the second is the total number of proposals
%                   for the drop proposals
%       ampMoves    a 1 x 2 vector where the first entry is the number of
%                   accepts and the second is the total number of proposals
%                   for the amplitude proposals
%       tau1Moves    a 1 x 2 vector where the first entry is the number of
%                   accepts and the second is the total number of proposals
%                   for the rise time proposals
%       tau2Moves    a 1 x 2 vector where the first entry is the number of
%                   accepts and the second is the total number of proposals
%                   for the fall time proposals

% number of samples in the trace
T = length(trace);

% unpack params struct for convenience
% see function get_params.m for information on the parameters for the
% algorithim

% params for priors
% rate
p_event=params.p_event;
% time constant bounds
tau1_min = params.tau1_min/params.dt;
tau1_max = params.tau1_max/params.dt;
tau2_min = params.tau2_min/params.dt;
tau2_max = params.tau2_max/params.dt;
% event amplitude bounds
a_min = params.a_min;
a_max = params.a_max;
% baseline bounds
b_min = params.b_min;
b_max = params.b_max;

% ar noise process params
% see Chib and Greenberg, Journal of Econometrics, 1994
% number of timesteps for filter
p = params.p;
% filter coefficient prior mean (MVN prior)
phi_0 = params.phi_0;
% precision/inverse covariance for filter coefficient prior (MVN)
Phi_0 = params.Phi_0; 
% params for inverse gamma prior on white-noise noise variance 
nu_0 = params.nu_0; 
sig2_0 = params.sig2_0; 
% if we don't want to use only a subset of the trace for noise estimation,
% then set subset to full trace
if isempty(params.noise_est_subset)
    params.noise_est_subset = 1:length(trace);
end

% parameters to balance sampling rates for different params
adddrop = params.add_drop_sweeps;
event_time_sweeps = params.event_time_sweeps;
amp_sweeps = params.amp_sweeps;
baseline_sweeps = params.baseline_sweeps;
tau1_sweeps = params.tau1_sweeps;
tau2_sweeps = params.tau2_sweeps;

% initial proposal variances for different params
time_proposal_var = params.time_proposal_var/params.dt;
tau1_std = params.tau1_prop_std/params.dt; %proposal variance of tau parameters
tau2_std = params.tau2_prop_std/params.dt; %proposal variance of tau parameters
a_std = params.amp_prop_std; %proposal variance of amplitude
b_std = params.baseline_prop_std; %propasal variance of baseline

% max length of event
event_samples = params.event_samples;
ef_init = genEfilt_ar([mean([tau1_min tau1_max]) mean([tau2_min tau2_max])],event_samples);


%dont let events get within x bins of eachother
exclusion_bound = params.exclusion_bound/params.dt;

% number of gibbs sweeps to do
num_sweeps = params.num_sweeps; 
burn_in_sweeps = params.burn_in_sweeps;
total_sweeps = num_sweeps + burn_in_sweeps;

% initialize arrays to hold samples and objective
samples_a  = [];
samples_b = [];
samples_s = [];
samples_tau_rise = [];
samples_tau_fall = [];
samples_phi = [];
samples_noise = [];
N_sto = [];
objective = [];

%initialize arrays for event features and noiseless trace
this_samp_amps = amps_init; % amplitudes
this_samp_times = times_init; % event times
this_samp_taus = cell(length(times_init),1); % array of lists of event this_samp_taus
phi = [1 zeros(1,p)];
efs = cell(1,times_init);
N = length(times_init);

% init some values
noiseless_trace = b*ones(1,T); %initial trace is set to baseline 
baseline = min(trace); %initial baseline value
N = length(this_samp_times); %number of events in trace
noise_var = params.noise_var_init;
residual = trace - noiseless_trace; %trace - prediction
% prior expected number of events for this trace
e_num_events = p_event*T;

% add initial events
for i = 1:length(times_init)
    
    % proposed tau and build exp filter
    this_samp_taus{i} = taus_init(i,:); 
    efs{i} = genEfilt_ar(this_samp_taus{i}, event_samples);
    
    % add event
    [~, noiseless_trace, residual] = ...
        add_event(times_init(1:i-1),noiseless_trace,residual,efs{i},...
        amps_init(i),this_samp_taus{i},trace,times_init(i),i);
    
end


% re-estimate the ar process parameters
if p > 0
    e = residual(params.noise_est_subset)'; % this is Tx1 (after transpose)
    E = [];
    for ip = 1:p
        E = [E e((p+1-ip):(end-ip))];
    end
    e = e((p+1):end);

    Phi_n = Phi_0 + noise_var^(-1)*(E'*E); %typo in paper
    phi_cond_mean = Phi_n\(Phi_0*phi_0 + noise_var^(-1)*E'*e);

    sample_phi = 1;
    while sample_phi
        phi = [1 mvnrnd(phi_cond_mean,inv(Phi_n))];

        phi_poly = -phi;
        phi_poly(1) = 1;
        if all(abs(roots(phi_poly))<1) %check stability
            sample_phi = 0;
        end
    end
end

% re-estimate the noise variance
df = (numel(noiseless_trace(params.noise_est_subset))); %DOF (possibly numel(ci(ti,:))-1)
d1 = -predAR(residual(params.noise_est_subset),phi,p,1 )/df; 
nu0 = nu_0; %nu_0 or 0
d0 = sig2_0; %sig2_0 or 0
A_samp = 0.5 * (df - p + nu0); %nu0 is prior
B_samp = 1/(0.5 * df * (d1 + d0)); %d0 is prior
noise_var = 1/gamrnd(A_samp,B_samp); %this could be inf but it shouldn't be


%% loop over sweeps to generate samples

addMoves = [0 0]; %first elem is number successful, second is number total
dropMoves = [0 0];
timeMoves = [0 0];
ampMoves = [0 0];
tau1Moves = [0 0];
tau2Moves = [0 0];

for i = 1:total_sweeps
    

    % do event time moves
    for ii = 1:event_time_sweeps
        
        for ni = 1:N%for each event
            proposed_time = this_samp_times(ni)+(time_proposal_var*randn); %add in noise 
            % bouncing off edges
            while proposed_time>T || proposed_time<0
                if proposed_time<0
                    proposed_time = -(proposed_time);
                elseif proposed_time>T
                    proposed_time = T-(proposed_time-T);
                end
            end
            
            %if its too close to another event, reject this move
            if any(abs(proposed_time-this_samp_times([1:(ni-1) (ni+1):end]))<exclusion_bound)
                continue
            end

            %create the proposal this_samp_times_tmp and noiseless_trace_tmp
            [this_samp_times_tmp, noiseless_trace_tmp, residual_tmp] = ...
                remove_event(this_samp_times,noiseless_trace,residual,...
                efs{ni},this_samp_amps(ni),this_samp_taus{ni},trace,this_samp_times(ni),ni);
            
            [this_samp_times_tmp, noiseless_trace_tmp, residual_tmp] = ...
                add_event(this_samp_times_tmp,noiseless_trace_tmp,...
                residual_tmp,efs{ni},this_samp_amps(ni),this_samp_taus{ni},trace,proposed_time,ni);

            %accept or reject
            prior_ratio = 1;
            ratio = exp(sum((1/(2*noise_var))*...
                ( predAR(residual_tmp,phi,p,1 ) - predAR(residual,phi,p,1 ) )))*prior_ratio;      
            
            if ratio>1 %accept
                this_samp_times = this_samp_times_tmp;
                noiseless_trace = noiseless_trace_tmp;
                residual = residual_tmp;
                timeMoves = timeMoves + [1 1];
                time_proposal_var = time_proposal_var + .1*rand*time_proposal_var/(i);
            elseif rand<ratio %accept
                this_samp_times = this_samp_times_tmp;
                noiseless_trace = noiseless_trace_tmp;
                residual = residual_tmp;
                timeMoves = timeMoves + [1 1];
                time_proposal_var = time_proposal_var + .1*rand*time_proposal_var/(i);
            else
                %reject - do nothing
                time_proposal_var = time_proposal_var - .1*rand*time_proposal_var/(i);
                timeMoves = timeMoves + [0 1];
            end
        end
    end

    
    % update amplitude of each event
    for ii = 1:amp_sweeps
        
        for ni = 1:N
            %sample with random walk proposal
            proposed_amp = this_samp_amps(ni)+(a_std*randn); %with bouncing off min and max
            while proposed_amp>a_max || proposed_amp<a_min
                if proposed_amp<a_min
                    proposed_amp = a_min+(a_min-proposed_amp);
                elseif proposed_amp>a_max
                    proposed_amp = a_max-(proposed_amp-a_max);
                end
            end

            % update sample with proposal
            [this_samp_times_tmp, noiseless_trace_tmp, residual_tmp] = ...
                remove_event(this_samp_times,noiseless_trace,residual,...
                efs{ni},this_samp_amps(ni),this_samp_taus{ni},trace,this_samp_times(ni),ni);
            [this_samp_times_tmp, noiseless_trace_tmp, residual_tmp] = ...
                add_event(this_samp_times_tmp,noiseless_trace_tmp,residual_tmp,...
                efs{ni},proposed_amp,this_samp_taus{ni},trace,this_samp_times(ni),ni);

            this_samp_amps_tmp = this_samp_amps;
            this_samp_amps_tmp(ni) = proposed_amp;

            %accept or reject - include a prior?
            prior_ratio = 1;
            ratio = exp(sum((1/(2*noise_var))*( predAR(residual_tmp,phi,p,1) - predAR(residual,phi,p,1) )))*prior_ratio;
            if ratio>1 %accept
                this_samp_amps = this_samp_amps_tmp;
                this_samp_times = this_samp_times_tmp;
                noiseless_trace = noiseless_trace_tmp;
                residual = residual_tmp;
                ampMoves = ampMoves + [1 1];
                a_std = a_std + 2*.1*rand*a_std/(i);
            elseif rand<ratio %accept
                this_samp_amps = this_samp_amps_tmp;
                this_samp_times = this_samp_times_tmp;
                noiseless_trace = noiseless_trace_tmp;
                residual = residual_tmp;
                ampMoves = ampMoves + [1 1];
                a_std = a_std + 2*.1*rand*a_std/(i);
            else
                %reject - do nothing
                a_std = a_std - .1*rand*a_std/(i);
                ampMoves = ampMoves + [0 1];
            end
        end
    end

    
    % update baseline of each trial
    for ii = 1:baseline_sweeps
        %sample with random walk proposal
        proposed_baseline = baseline+(b_std*randn); %with bouncing off min and max
        while proposed_baseline>b_max || proposed_baseline<b_min
            if proposed_baseline<b_min
                proposed_baseline = b_min+(b_min-proposed_baseline);
            elseif proposed_baseline>b_max
                proposed_baseline = b_max-(proposed_baseline-b_max);
            end
        end

        % update with proposal
        [noiseless_trace_tmp, residual_tmp] = ...
            remove_base(noiseless_trace,residual,baseline,trace);   
        [noiseless_trace_tmp, residual_tmp] = ...
            add_base(noiseless_trace_tmp,residual_tmp,proposed_baseline,trace);

        %accept or reject - include a prior?
        prior_ratio = 1;
        ratio = exp(sum((1/(2*noise_var))*(predAR(residual_tmp,phi,p,1 ) - predAR(residual,phi,p,1 ))))*prior_ratio;
        if ratio>1 %accept
            baseline = proposed_baseline;
            noiseless_trace = noiseless_trace_tmp;
            residual = residual_tmp;
            b_std = b_std + 2*.1*rand*b_std/(i);
        elseif rand<ratio %accept
            baseline = proposed_baseline;
            noiseless_trace = noiseless_trace_tmp;
            residual = residual_tmp;
            b_std = b_std + 2*.1*rand*b_std/(i);
        else
            %reject - do nothing
            b_std = b_std - .1*rand*b_std/(i);
        end
    end

    if i>1
        %% this is the section that updates the number of events (add/drop)
        % loop over add/drop a few times
        % define insertion proposal distribution as the likelihood function
        % define removal proposal distribution as uniform over events
        % perhaps better is to choose smarter removals.
        for ii = 1:adddrop 
            %propose a uniform add
            %pick a random point
            proposed_time = T*rand;
            %dont add if we have too many events or the proposed new location
            %is too close to another one
            if ~any(abs(proposed_time-this_samp_times)<exclusion_bound)
                
                %must add event to each trial (at mean location or sampled -- more appropriate if sampled, but make sure no trial's event violates exclusion)
                residual_tmp = residual;
                
                % initialize amplitude based on value of trace around event
                start_ind = max(1,floor(proposed_time));
                end_ind = min(start_ind + params.a_init_window*2,length(residual));
                [local_max,tmpi_tmp] = max(residual(start_ind:end_ind));
                proposed_time = tmpi_tmp + start_ind - 1;
                a_init = max(local_max + a_std*randn,a_min);
                
                [this_samp_times_tmp, noiseless_trace_tmp, residual_tmp] =...
                    add_event(this_samp_times,noiseless_trace,residual_tmp,...
                    ef_init,a_init,tau,trace,proposed_time, N+1);
                this_samp_amps_tmp = [this_samp_amps a_init];
                
                %accept or reject
                % prior probs
                fprob = 1;%1/T(1);%forward probability
                rprob = 1;%1/(N+1);%reverse (remove at that spot) probability
                ratio = exp(sum((1./(2*noise_var)).*...
                    ( predAR(residual_tmp,phi,p,1 ) - predAR(residual,phi,p,1 ) )))*...
                    (rprob/fprob)*e_num_events/(N+1);
                if (ratio>1)||(ratio>rand) %accept
                    this_samp_amps = this_samp_amps_tmp;
                    this_samp_times = this_samp_times_tmp;
                    noiseless_trace = noiseless_trace_tmp;
                    this_samp_taus{N+1} = tau;
                    efs{N+1} = ef_init;
                    residual = residual_tmp;
                    addMoves = addMoves + [1 1];
                else
                    %reject - do nothing
                    addMoves = addMoves + [0 1];
                end
                N = length(this_samp_times);
            end


            % delete
            if N>0%i.e. we if have at least one event 
                
                %propose a uniform removal
                proposed_event_i = randi(N);%pick one of the events at random
                %must remove event from each trial
                residual_tmp = residual;
                %always remove the ith event (the ith event of each trial is linked)                     
                [this_samp_times_tmp, noiseless_trace_tmp, residual_tmp] =...
                    remove_event(this_samp_times,noiseless_trace,residual_tmp,...
                    efs{proposed_event_i},this_samp_amps(proposed_event_i),...
                    this_samp_taus{proposed_event_i},trace,this_samp_times(proposed_event_i),...
                    proposed_event_i);                

                %accept or reject
                %reverse probability
                rprob = 1;%1/T(1);
                %compute forward prob
                fprob = 1;%1/N;
                %posterior times reverse prob/forward prob
                ratio = exp(sum((1./(2*noise_var)).*( predAR(residual_tmp,phi,p,1 ) - predAR(residual,phi,p,1 ) )))*(rprob/fprob)*N/e_num_events;%((T(1)-e_num_events(1))/e_num_events(1)); 
                if (ratio>1)||(ratio>rand)%accept
                    this_samp_amps(proposed_event_i) = [];
                    this_samp_times = this_samp_times_tmp;
                    noiseless_trace = noiseless_trace_tmp;
                    this_samp_taus(proposed_event_i) = [];
                    efs(proposed_event_i) = [];
                    residual = residual_tmp;
                    dropMoves = dropMoves + [1 1]; 
                else
                    %reject - do nothing
                    dropMoves = dropMoves + [0 1];
                end
                N = length(this_samp_times);
            end
        end
    end
    
 
    %% this is the section that updates rise tau
    % update tau (via random walk sampling)   
    for ii = 1:tau1_sweeps
        for ni = 1:N 
            
            % propose new rise tau
            proposed_tau = this_samp_taus{ni};
            proposed_tau(1) = proposed_tau(1)+(tau1_std*randn); %with bouncing off min and max
            tau_max = min([proposed_tau(2) tau1_max]);
            tau_min = tau1_min;
            while proposed_tau(1)>tau_max || proposed_tau(1)<tau_min
                if proposed_tau(1) < tau_min
                    proposed_tau(1) = tau_min+(tau_min-proposed_tau(1));
                elseif proposed_tau(1)>tau_max
                    proposed_tau(1) = tau_max -(proposed_tau(1)-tau_max);
                end
            end 

            ef_ = genEfilt_ar(proposed_tau,event_samples);%exponential filter

            % udpate with proposal   
            [~, noiseless_trace_tmp, residual_tmp] = ...
                remove_event(this_samp_times,noiseless_trace,residual,...
                efs{ni},this_samp_amps(ni),this_samp_taus{ni},trace,...
                this_samp_times(ni),ni);
            [~, noiseless_trace_tmp, residual_tmp] = ...
                add_event(this_samp_times,noiseless_trace_tmp,residual_tmp,...
                ef_,this_samp_amps(ni),proposed_tau,trace,this_samp_times(ni),ni);

            %accept or reject
            prior_ratio = 1;
%             prior_ratio = (gampdf(proposed_tau(1),1.5,1))/(gampdf(tau(1),1.5,1));
            ratio = exp(sum(sum((1./(2*noise_var)).*...
                ( predAR(residual_tmp,phi,p,1 ) - predAR(residual,phi,p,1 ) ))))*prior_ratio;
            if ratio>1 %accept
                noiseless_trace = noiseless_trace_tmp;
                residual = residual_tmp;
                this_samp_taus{ni} = proposed_tau;
                efs{ni} = ef_;
                tau1Moves = tau1Moves + [1 1];
                tau1_std = tau1_std + .1*rand*tau1_std/(i);
            elseif rand<ratio %accept
                noiseless_trace = noiseless_trace_tmp;
                residual = residual_tmp;
                this_samp_taus{ni} = proposed_tau;
                efs{ni} = ef_;
                tau1Moves = tau1Moves + [1 1];
                tau1_std = tau1_std + .1*rand*tau1_std/(i);
            else
                %reject - do nothing
                tau1_std = tau1_std - .1*rand*tau1_std/(i);
                tau1Moves = tau1Moves + [0 1];
            end
        end
    end
    

%% this is the section that updates fall tau

    % update tau (via random walk sampling)
    for ii = 1:tau2_sweeps
        for ni = 1:N 
            % update both tau values
            proposed_tau = this_samp_taus{ni};    
            proposed_tau(2) = proposed_tau(2)+(tau2_std*randn);
            tau_min = max([proposed_tau(1) tau2_min]);
            tau_max = tau2_max;
            while proposed_tau(2)>tau_max || proposed_tau(2)<proposed_tau(1)
                if proposed_tau(2)<tau_min
                    proposed_tau(2) = tau_min+(tau_min-proposed_tau(2));
                elseif proposed_tau(2)>tau_max
                    proposed_tau(2) = tau_max-(proposed_tau(2)-tau_max);
                end
            end  
            
            ef_ = genEfilt_ar(proposed_tau,event_samples);%exponential filter

            % update proposal
            [~, noiseless_trace_tmp, residual_tmp] = ...
                remove_event(this_samp_times,noiseless_trace,residual,...
                efs{ni},this_samp_amps(ni),this_samp_taus{ni},trace,this_samp_times(ni),ni);
            [~, noiseless_trace_tmp, residual_tmp] = ...
                add_event(this_samp_times,noiseless_trace_tmp,residual_tmp,...
                ef_,this_samp_amps(ni),proposed_tau,trace,this_samp_times(ni),ni);

            %accept or reject
            prior_ratio = 1;
%             prior_ratio = gampdf(proposed_tau(2),12,1)/gampdf(tau(2),12,1);
            ratio = exp(sum(sum((1./(2*noise_var)).*( predAR(residual_tmp,phi,p,1 ) - predAR(residual,phi,p,1 ) ))))*prior_ratio;
            if ratio>1 %accept
                noiseless_trace = noiseless_trace_tmp;
                residual = residual_tmp;
                this_samp_taus{ni} = proposed_tau;
                efs{ni} = ef_;
                tau2Moves = tau2Moves + [1 1];
                tau2_std = tau2_std + .1*rand*tau2_std/(i);
            elseif rand<ratio %accept
                noiseless_trace = noiseless_trace_tmp;
                residual = residual_tmp;
                this_samp_taus{ni} = proposed_tau;
                efs{ni} = ef_;
                tau2Moves = tau2Moves + [1 1];
                tau2_std = tau2_std + .1*rand*tau2_std/(i);
            else
                %reject - do nothing
                tau2_std = tau2_std - .1*rand*tau2_std/(i);
                tau2Moves = tau2Moves + [0 1];
            end
        end
    end

    
    % re-estimate the ar process parameters
    if p>0 %&& i>(num_sweeps/100)
        e = residual(params.noise_est_subset)'; % this is Tx1 (after transpose)
        E = [];
        for ip = 1:p
            E = [E e((p+1-ip):(end-ip))];
        end
        e = e((p+1):end);

        Phi_n = Phi_0 + noise_var^(-1)*(E'*E); %typo in paper

        phi_cond_mean = Phi_n\(Phi_0*phi_0 + noise_var^(-1)*E'*e);

%         keyboard
        sample_phi = 1;
        while sample_phi
            phi = [1 mvnrnd(phi_cond_mean,inv(Phi_n))];

            phi_poly = -phi;
            phi_poly(1) = 1;
            if all(abs(roots(phi_poly))<1) %check stability
                sample_phi = 0;
            end
        end
        
    end

    % re-estimate the noise variance
    df = (numel(noiseless_trace(params.noise_est_subset))); %DOF (possibly numel(ci(ti,:))-1)
    d1 = -predAR(residual(params.noise_est_subset),phi,p,1 )/df; 
    nu0 = nu_0; %nu_0 or 0
    d0 = sig2_0; %sig2_0 or 0

    A_samp = 0.5 * (df - p + nu0); %nu0 is prior
    B_samp = 1/(0.5 * df * (d1 + d0)); %d0 is prior
    noise_var = 1/gamrnd(A_samp,B_samp); %this could be inf but it shouldn't be


    %store samples
    N_sto = [N_sto N];
    samples_a = [samples_a this_samp_amps]; %trial amplitudes
    samples_b = [samples_b baseline]; %trial baselines
    samples_s = [samples_s this_samp_times]; %shared events
    for event_i = 1:N
        samples_tau_rise = [samples_tau_rise this_samp_taus{event_i}(1)]; %save tau values
        samples_tau_fall = [samples_tau_fall this_samp_taus{event_i}(2)];
    end
    samples_phi = [samples_phi; phi];
    samples_noise = [samples_noise noise_var];
    
    %store overall logliklihood as well
    objective = [objective -T/2*log(noise_var) + predAR(residual,phi,p,1 )/(2*noise_var) + N*log(e_num_events) - log(factorial(N))];

end

%% build output structs

%details about what the mcmc did
%addMoves, dropMoves, and timeMoves give acceptance probabilities for each subclass of move
mcmc.addMoves=addMoves;
mcmc.timeMoves=timeMoves;
mcmc.dropMoves=dropMoves;
mcmc.ampMoves=ampMoves;
mcmc.tau1Moves=tau1Moves;
mcmc.tau2Moves=tau2Moves;


posterior.amp=samples_a;
posterior.base=samples_b;
posterior.tau1=samples_tau_rise;
posterior.tau2=samples_tau_fall;
posterior.num_events = N_sto;
posterior.phi=samples_phi;
posterior.noise = samples_noise;
posterior.obj = objective;
posterior.times = samples_s;






 

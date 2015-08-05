function [trials, mcmc, params]  = sampleParams_ARnoise_splittau(trace,tau, Tguess, params)
%parameters


%noise level here matters for the proposal distribution (how much it 
%should trust data for proposals vs how much it should mix based on uniform prior)
%this is accounted for by calciumNoiseVar
NoiseVar_init=5; %inial noise estimate
% p_spike=1/40;%what percent of the bins hacve a spike in then
p_spike=params.p_spike;
proposalVar=10;
nsweeps=1000; %number of sweeps of sampler
% if acceptance rates are too high, increase proposal width, 
% if too low, decrease them (for time moves, tau, amplitude)
% tau_std = 1;
tau1_std = 2/20000/params.dt; %proposal variance of tau parameters
tau2_std = 2/20000/params.dt; %proposal variance of tau parameters
tau1_min = params.tau1_min;
tau1_max = params.tau1_max;
tau2_min = params.tau2_min;
tau2_max = params.tau2_max;
%all of these are multiplied by big A
a_std = .2; %proposal variance of amplitude
a_min = params.a_min;
a_max = Inf;
b_std = 2; %propasal variance of baseline  % was at 0.3 ... increased it for in vivo data
b_min = -50;
b_max = 50;
exclusion_bound = 1;%dont let bursts get within x bins of eachother. this should be in time
% maxNbursts = 3;%if we want to add bursts, whats the maximum bnumber that we will look for?

Dt=1; %bin unit - don't change this
% A=400; % scale factor for all magnitudes for this calcium data setup
A=1; % scale factor for all magnitudes for this calcium data setup
% b=0; %initial baseline value
b=min(trace); %initial baseline value
nu_0 = 5; %prior on shared burst time - ntrials
sig2_0 = .1; %prior on shared burst time - variance

%noise model is AR(p)
if isfield(params, 'p')
    p = params.p;
else
    p = 4;
end
% phi prior
phi_0 = zeros(p,1);
Phi_0 = 10*eye(p); %inverse covariance 3

adddrop = 5;
% maxNbursts = length(Tguess);
maxNbursts = Inf;

indreport=.1:.1:1;
indreporti=round(nsweeps*indreport);
fprintf('Progress:')

% initialize some parameters
nBins = length(trace); %for all of this, units are bins and spiketrains go from 0 to T where T is number of bins
fBins = 2000;
ef = genEfilt_ar(tau,fBins);%exponential filter
ef_init = ef;

samples_a  = cell(1,nsweeps);
samples_b = cell(1,nsweeps);
samples_s = cell(1,nsweeps);
samples_pr = cell(1,nsweeps);
samples_tau = cell(1,nsweeps);
samples_phi = cell(1,nsweeps);
samples_noise = cell(1,nsweeps);
N_sto = [];
objective = [];

NoiseVar = NoiseVar_init; %separate calcium per trial
baseline = b;

% intiailize burst train and predicted calcium
%this is based on simply what we tell it. 

%initialize spikes and calcium
ati = []; % array of lists of spike times
sti = []; % array of lists of spike times
sti_ = []; % array of lists of spike times
taus = cell(1); % array of lists of event taus
phi = [1 zeros(1,p)];

efs = cell(1,length(Tguess));

pr = b*ones(1,nBins); %initial calcium is set to baseline 

N = length(sti); %number of spikes in spiketrain

%initial logC - compute likelihood initially completely - updates to likelihood will be local
%for AR(p) noise, we need a different difference inside the likelihood
diffY = (trace-pr); %trace - prediction

m = p_spike*nBins;

diffY_ = diffY;
for i = 1:length(Tguess)
    efs{i} = ef;
    tmpi = Tguess(i); 
    sti_ = [sti tmpi];
    %must add spike to each trial (at mean location or sampled -- more appropriate if sampled)
    pr_ = pr;
    ati_ = ati;
    a_init = max(trace(tmpi)/A,a_min);
    [sti_, pr_, diffY_] = addSpike_ar(sti,pr,diffY_,efs{i},a_init,tau,trace,tmpi, N+1, Dt, A); %adds all trials' spikes at same time
    taus{i} = tau;
    ati_ = [ati_ a_init];
    ati = ati_;
    sti = sti_;
    pr = pr_;
    N = length(sti); %number of spikes in spiketrain
end
diffY = diffY_;

sti_= sti;
diffY_= diffY;
N=length(sti);

%% loop over sweeps to generate samples
addMoves = [0 0]; %first elem is number successful, second is number total
dropMoves = [0 0];
timeMoves = [0 0];
ampMoves = [0 0];
tauMoves = [0 0];
for i = 1:nsweeps
    
%     if mod(i,10) == 0
%         disp(length(ati))
%     end
    
    % do burst time moves
    for ii = 1:3
        %guess on time and amplitude
        si = sti;
        ai = ati;
        for ni = 1:N%for each burst
            tmpi = si(ni);
            tmpi_ = si(ni)+(proposalVar*randn); %add in noise 
            % bouncing off edges
            while tmpi_>nBins || tmpi_<0
                if tmpi_<0
                    tmpi_ = -(tmpi_);
                elseif tmpi_>nBins
                    tmpi_ = nBins-(tmpi_-nBins);
                end
            end
            %if its too close to another burst, reject this move
            if any(abs(tmpi_-si([1:(ni-1) (ni+1):end]))<exclusion_bound)
                continue
            end

            %create the proposal si_ and pr_
            %update logC_ to adjusted
            [si_, pr_, diffY_] = removeSpike_ar(si,pr,diffY,efs{ni},ai(ni),taus{ni},trace,tmpi,ni, Dt, A);
            [si_, pr_, diffY_] = addSpike_ar(si_,pr_,diffY_,efs{ni},ai(ni),taus{ni},trace,tmpi_,ni, Dt, A);

            %accept or reject
            %for prior: (1) use ratio or(2) set prior to 1.
            prior_ratio = 1;

            ratio = exp(sum((1/(2*NoiseVar))*( predAR(diffY_,phi,p,1) - predAR(diffY,phi,p,1) )))*prior_ratio;            
            if ratio>1 %accept
                si = si_;
                pr = pr_;
                diffY = diffY_;
                timeMoves = timeMoves + [1 1];
                proposalVar = proposalVar + .1*rand*proposalVar/(i);
            elseif rand<ratio %accept
                si = si_;
                pr = pr_;
                diffY = diffY_;
                timeMoves = timeMoves + [1 1];
                proposalVar = proposalVar + .1*rand*proposalVar/(i);
            else
                %reject - do nothing
                proposalVar = proposalVar - .1*rand*proposalVar/(i);
                timeMoves = timeMoves + [0 1];
            end
        end
        sti = si;
    end

    
    % update amplitude of each burst
    for ii = 1:5
        si = sti; 
        ai = ati;
        for ni = 1:N
            %sample with random walk proposal
            tmp_a = ai(ni);
            tmp_a_ = tmp_a+(a_std*randn); %with bouncing off min and max
            while tmp_a_>a_max || tmp_a_<a_min
                if tmp_a_<a_min
                    tmp_a_ = a_min+(a_min-tmp_a_);
                elseif tmp_a_>a_max
                    tmp_a_ = a_max-(tmp_a_-a_max);
                end
            end

            %set si_ to set of bursts with the move and pr_ to adjusted calcium and update logC_ to adjusted
            [si_, pr_, diffY_] = removeSpike_ar(si,pr,diffY,efs{ni},ai(ni),taus{ni},trace,si(ni),ni, Dt, A);
            [si_, pr_, diffY_] = addSpike_ar(si_,pr_,diffY_,efs{ni},tmp_a_,taus{ni},trace,si(ni),ni, Dt, A);

            ai_ = ai;
            ai_(ni) = tmp_a_;

            %accept or reject - include a prior?
            prior_ratio = 1;
            ratio = exp(sum((1/(2*NoiseVar))*( predAR(diffY_,phi,p,1) - predAR(diffY,phi,p,1) )))*prior_ratio;
            if ratio>1 %accept
                ai = ai_;
                si = si_;
                pr = pr_;
                diffY = diffY_;
                ampMoves = ampMoves + [1 1];
                a_std = a_std + 2*.1*rand*a_std/(i);
            elseif rand<ratio %accept
                ai = ai_;
                si = si_;
                pr = pr_;
                diffY = diffY_;
                ampMoves = ampMoves + [1 1];
                a_std = a_std + 2*.1*rand*a_std/(i);
            else
                %reject - do nothing
                a_std = a_std - .1*rand*a_std/(i);
                ampMoves = ampMoves + [0 1];
            end
        end
        ati = ai;
    end

    
    % update baseline of each trial
    for ii = 1:1
        %sample with random walk proposal
        tmp_b = baseline;
        tmp_b_ = tmp_b+(b_std*randn); %with bouncing off min and max
        while tmp_b_>b_max || tmp_b_<b_min
            if tmp_b_<b_min
                tmp_b_ = b_min+(b_min-tmp_b_);
            elseif tmp_b_>b_max
                tmp_b_ = b_max-(tmp_b_-b_max);
            end
        end

        %set si_ to set of bursts with the move and pr_ to adjusted calcium and update logC_ to adjusted
        [pr_, diffY_] = remove_base_ar(pr,diffY,tmp_b,trace,A);   
        [pr_, diffY_] = add_base_ar(pr_,diffY_,tmp_b_,trace,A);

        %accept or reject - include a prior?
        prior_ratio = 1;
        ratio = exp(sum((1/(2*NoiseVar))*(  predAR(diffY_,phi,p,1) - predAR(diffY,phi,p,1)  )))*prior_ratio;
        if ratio>1 %accept
            baseline = tmp_b_;
            pr = pr_;
            diffY = diffY_;
            b_std = b_std + 2*.1*rand*b_std/(i);
        elseif rand<ratio %accept
            baseline = tmp_b_;
            pr = pr_;
            diffY = diffY_;
            b_std = b_std + 2*.1*rand*b_std/(i);
        else
            b_std = b_std - .1*rand*b_std/(i);
            %reject - do nothing
        end
    end

    if i>1
    %% this is the section that updates the number of spikes (add/drop)
    % loop over add/drop a few times
    %define insertion proposal distribution as the likelihood function
    %define removal proposal distribution as uniform over bursts
    %perhaps better is to choose smarter removals.
        for ii = 1:adddrop 
            %propose a uniform add
            %pick a random point
            tmpi = min(nBins)*rand;
            %dont add if we have too many bursts or the proposed new location
            %is too close to another one
            if ~(any(abs(tmpi-sti)<exclusion_bound) || N >= maxNbursts)
                sti_ = [sti tmpi];
                %must add burst to each trial (at mean location or sampled -- more appropriate if sampled, but make sure no trial's burst violates exclusion)
                diffY_ = diffY;
                pr_ = pr;
                ati_ = ati;
                a_init = max(trace(max(1,floor(tmpi)))/A - baseline + a_std*randn,a_min);%propose an initial amplitude for it
                [si_, pr_, diffY_] = addSpike_ar(sti,pr,diffY_,ef_init,a_init,tau,trace,tmpi, N+1, Dt, A); %adds all trials' bursts at same time
                sti_ = si_;
                ati_ = [ati_ a_init];
                fprob = 1/nBins(1);%forward probability
                rprob = 1/(N+1);%reverse (remove at that spot) probability
                %accept or reject
%                 figure(120)
%                 plot(trace)
%                 hold on
%                 plot(pr_,'r')
%                 hold off
%                 drawnow
%                 pause
                ratio = exp(sum((1./(2*NoiseVar)).*( predAR(diffY_,phi,p,1) - predAR(diffY,phi,p,1) )))*(rprob/fprob)*(m(1)/(nBins(1)-m(1))); %posterior times reverse prob/forward prob
                if (ratio>1)||(ratio>rand) %accept
                    ati = ati_;
                    sti = sti_;
                    pr = pr_;
                    taus{N+1} = tau;
                    efs{N+1} = ef_init;
                    diffY = diffY_;
                    addMoves = addMoves + [1 1];
                else
                    %reject - do nothing
                    addMoves = addMoves + [0 1];
                end
                N = length(sti);
            end


            % delete
            if N>0%i.e. we if have at least one spike           
                %propose a uniform removal
                tmpi = randi(N);%pick one of the spikes at random
                sti_ = sti;
                sti_(tmpi) = [];
                %must remove burst from each trial
                diffY_ = diffY;
                pr_ = pr;
                ati_ = ati;
                %always remove the ith burst (the ith burst of each trial is linked)                     
                [si_, pr_, diffY_] = removeSpike_ar(sti,pr,diffY_,efs{tmpi},ati(tmpi),taus{tmpi},trace,sti(tmpi),tmpi, Dt, A);
                sti_ = si_;
                ati_(tmpi) = [];

                %reverse probability
                rprob = 1/nBins(1);

                %compute forward prob
                fprob = 1/N;

                %accept or reject
                %posterior times reverse prob/forward prob
                ratio = exp(sum((1./(2*NoiseVar)).*( predAR(diffY_,phi,p,1) - predAR(diffY,phi,p,1) )))*(rprob/fprob)*((nBins(1)-m(1))/m(1)); 
                if (ratio>1)||(ratio>rand)%accept
                    ati = ati_;
                    sti = sti_;
                    pr = pr_;
                    taus(tmpi) = [];
                    efs(tmpi) = [];
                    diffY = diffY_;
                    dropMoves = dropMoves + [1 1]; 
                else
                    %reject - do nothing
                    dropMoves = dropMoves + [0 1];
                end
                N = length(sti);
            end
        end
    end
    
    
    for ii = 1:1
        for ni = 1:N 
            % update both tau values
            tau_ = taus{ni};
            tau_(1) = tau_(1)+(tau1_std*randn); %with bouncing off min and max
            tau_max = min([tau(2) tau1_max]);
            tau_min = tau1_min;
            while tau_(1)>tau_max || tau_(1)<tau_min
                if tau_(1) < tau_min
                    tau_(1) = tau_min+(tau_min-tau_(1));
                elseif tau_(1)>tau_max
                    tau_(1) = tau_max -(tau_(1)-tau_max);
                end
            end 

            ef_ = genEfilt_ar(tau_,fBins);%exponential filter

            %remove all old bumps and replace them with new bumps    
            diffY_ = diffY;
            pr_ = pr;
            [~, pr_, diffY_] = removeSpike_ar(sti,pr_,diffY_,efs{ni},ati(ni),taus{ni},trace,sti(ni),ni, Dt, A);
            [~, pr_, diffY_] = addSpike_ar(sti,pr_,diffY_,ef_,ati(ni),tau_,trace,sti(ni),ni, Dt, A);

            %accept or reject
            prior_ratio = 1;
%             prior_ratio = (gampdf(tau_(1),1.5,1))/(gampdf(tau(1),1.5,1));
            ratio = exp(sum(sum((1./(2*NoiseVar)).*( predAR(diffY_,phi,p,1) - predAR(diffY,phi,p,1) ))))*prior_ratio;
            if ratio>1 %accept
                pr = pr_;
                diffY = diffY_;
                taus{ni} = tau_;
                efs{ni} = ef_;
                tauMoves = tauMoves + [1 1];
               tau1_std = tau1_std + .1*rand*tau1_std/(i);
            elseif rand<ratio %accept
                pr = pr_;
                diffY = diffY_;
                taus{ni} = tau_;
                efs{ni} = ef_;
                tauMoves = tauMoves + [1 1];
               tau1_std = tau1_std + .1*rand*tau1_std/(i);
            else
                %reject - do nothing
               tau1_std = tau1_std - .1*rand*tau1_std/(i);
                tauMoves = tauMoves + [0 1];
            end
        end
    end
    

     %% this is the section that updates tau
    % update tau (via random walk sampling)
    for ii = 1:1  
        for ni = 1:N 
            % update both tau values
            tau_ = taus{ni};    
            tau_(2) = tau_(2)+(tau2_std*randn);
            tau_min = max([tau(1) tau2_min]);
            tau_max = tau2_max;
            while tau_(2)>tau_max || tau_(2)<tau_(1)
                if tau_(2)<tau_min
                    tau_(2) = tau_min+(tau_min-tau_(2));
                elseif tau_(2)>tau_max
                    tau_(2) = tau_max-(tau_(2)-tau_max);
                end
            end  

            ef_ = genEfilt_ar(tau_,fBins);%exponential filter

            %remove all old bumps and replace them with new bumps    
            diffY_ = diffY;
            pr_ = pr;
            [~, pr_, diffY_] = removeSpike_ar(sti,pr_,diffY_,efs{ni},ati(ni),taus{ni},trace,sti(ni),ni, Dt, A);
            [~, pr_, diffY_] = addSpike_ar(sti,pr_,diffY_,ef_,ati(ni),tau_,trace,sti(ni),ni, Dt, A);

            %accept or reject
            prior_ratio = 1;
%             prior_ratio = gampdf(tau_(2),12,1)/gampdf(tau(2),12,1);
            ratio = exp(sum(sum((1./(2*NoiseVar)).*( predAR(diffY_,phi,p,1) - predAR(diffY,phi,p,1) ))))*prior_ratio;
            if ratio>1 %accept
                pr = pr_;
                diffY = diffY_;
                taus{ni} = tau_;
                efs{ni} = ef_;
                tauMoves = tauMoves + [1 1];
                tau2_std = tau2_std + .1*rand*tau2_std/(i);
            elseif rand<ratio %accept
                pr = pr_;
                diffY = diffY_;
                taus{ni} = tau_;
                efs{ni} = ef_;
                tauMoves = tauMoves + [1 1];
                tau2_std = tau2_std + .1*rand*tau2_std/(i);
            else
                %reject - do nothing
                tau2_std = tau2_std - .1*rand*tau2_std/(i);
                tauMoves = tauMoves + [0 1];
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % re-estimate the noise model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % re-estimate the process parameters
    %%%%%%%%%%%%%%%%
    % estimate phi (ignore initial condition boundary effects)
    %%%%%%%%%%%%%%%%
    if p>0 %&& i>(nsweeps/100)
        e = diffY'; % this is Tx1 (after transpose)
        E = [];
        for ip = 1:p
            E = [E e((p+1-ip):(end-ip))];
        end
        e = e((p+1):end);

        Phi_n = Phi_0 + NoiseVar^(-1)*(E'*E); %typo in paper

        phi_cond_mean = Phi_n\(Phi_0*phi_0 + NoiseVar^(-1)*E'*e);

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
    
    %%%%%%%%%%%%%%%%%%%%%
    % estimate noise
    %%%%%%%%%%%%%%%%%%%%%
    % re-estimate the noise variance
    if ~isempty(sti)
        df = (numel(pr)); %DOF (possibly numel(ci(ti,:))-1)
        d1 = -predAR(diffY,phi,p,1)/df; 
        nu0 = nu_0; %nu_0 or 0
        d0 = sig2_0; %sig2_0 or 0
        
        A_samp = 0.5 * (df + nu0); %nu0 is prior
        B_samp = 1/(0.5 * df * (d1 + d0)); %d0 is prior
        NoiseVar = 1/gamrnd(A_samp,B_samp); %this could be inf but it shouldn't be
%         if ~isfinite(NoiseVar)
%             keyboard
%         end
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    %store things
    N_sto = [N_sto N];
    samples_a{i} = ati; %trial amplitudes
    samples_b{i} = baseline; %trial baselines
    samples_s{i} = sti; %shared bursts
    samples_pr{i} = pr; %save calcium traces
    samples_tau{i} = taus; %save tau values
    samples_phi{i} = phi;
    samples_noise{i} = NoiseVar;
    %store overall logliklihood as well
%     if abs(sum(logC)-sum(sum(-(pr)-cell2mat(trace)).^2)))>1
%         figure(90)
%         subplot(121)
%         plot(cell2mat(samples_c{i-1})')
%         subplot(122)
%         plot(cell2mat(pr)')
%         keyboard
%     end

    objective = [objective -predAR(diffY,phi,p,1)];
%     figure(10);
%     plot(diffY_)
%     drawnow
%     plot(ci{1});hold on;
%     plot(CaF{1},'r');hold off
    if sum(ismember(indreporti,i))
        fprintf([num2str(indreport(ismember(indreporti,i)),2),', '])
    end
end

%% Vigi's Clean up
%details about what the mcmc did
%addMoves, dropMoves, and timeMoves give acceptance probabilities for each subclass of move
mcmc.addMoves=addMoves;
mcmc.timeMoves=timeMoves;
mcmc.dropMoves=dropMoves;
mcmc.ampMoves=ampMoves;
mcmc.tauMoves=tauMoves;
mcmc.N_sto=N_sto;%number of bursts

trials.amp=samples_a;
trials.base=samples_b;
% trials.curves=samples_pr;
trials.tau=samples_tau;
trials.phi=samples_phi;
trials.noise = samples_noise;
trials.obj = objective;
trials.times = samples_s;

params.NoiseVar = NoiseVar_init; 
params.proposalVar = proposalVar;
params.nsweeps = nsweeps;
params.tau1_std = tau1_std; %proposal variance of tau parameters
params.tau2_std = tau2_std; %proposal variance of tau parameters
params.tau1_min = tau1_min;
params.tau1_max = tau1_max;
params.tau2_min = tau2_min;
params.tau2_max = tau2_max;
params.tau_min = tau1_min;
params.tau_max = tau2_max;
params.a_std = a_std; %proposal variance of amplitude
params.a_min = a_min;
params.a_max = a_max;
params.b_std = b_std; %propasal variance of baseline
params.b_min = b_min;
params.b_max = b_max;
params.exclusion_bound = exclusion_bound;
params.Dt=Dt; %bin unit - don't change this
params.A=A; % scale factor for all magnitudes for this calcium data setup
params.b=b; %initial baseline value
params.p_spike = p_spike;
params.p = p;
params.phi_0 = phi_0;
params.Phi_0 = Phi_0;

% disp('Below are the moves that were done')
% display(['time: ' num2str(timeMoves(1)/timeMoves(2))])
% display(['add: ' num2str(addMoves(1)/addMoves(2))])
% display(['drop: ' num2str(dropMoves(1)/dropMoves(2))])
% display(['amplitude: ' num2str(ampMoves(1)/ampMoves(2))])
% display(['tau: ' num2str(tauMoves(1)/tauMoves(2))])




function [trials, mcmc, params]  = sampleParams(trace,tau, Tguess,params)
%parameters
%noise level here matters for the proposal distribution (how much it 
%should trust data for proposals vs how much it should mix based on uniform prior)
%this is accounted for by calciumNoiseVar
NoiseVar_init=5; %inial noise estimate
% p_spike=1/40;%what percent of the bins hacve a spike in then
% p_spike=1/100000000000000000000000000000000;%what percent of the bins hacve a spike in then
p_spike=params.p_spike;
proposalVar=100000;%likeliness to accept moves
nsweeps=200; %number of sweeps of sampler
% if acceptance rates are too high, increase proposal width, 
% if too low, decrease them (for time moves, tau, amplitude)
% tau_std = 1;
tau1_std = 2/20000/params.dt; %proposal variance of tau parameters
tau2_std = 2/20000/params.dt; %proposal variance of tau parameters
tau_min = params.tau_min;%params.tau1_min;
tau_max = params.tau_max;%params.tau1_max;
%all of these are multiplied by big A
a_std = .2; %proposal variance of amplitude
a_min = params.a_min;
a_max = Inf;
b_std = .3; %propasal variance of baseline
b_min = 0;
b_max = 30;
exclusion_bound = 1;%dont let bursts get within x bins of eachother. this should be in time
% maxNbursts = 3;%if we want to add bursts, whats the maximum bnumber that we will look for?

Dt=1; %bin unit - don't change this
% A=400; % scale factor for all magnitudes for this calcium data setup
A=1; % scale factor for all magnitudes for this calcium data setup
b=median(trace); %initial baseline value
nu_0 = 5; %prior on shared burst time - ntrials
sig2_0 = .1; %prior on shared burst time - variance

adddrop = 5;
% maxNbursts = length(Tguess);
maxNbursts = length(Tguess) + 8;
maxNbursts = Inf;

indreport=.1:.1:1;
indreporti=round(nsweeps*indreport);
% fprintf('Progress:')


% initialize some parameters
nBins = length(trace); %for all of this, units are bins and spiketrains go from 0 to T where T is number of bins
fBins = 2000;
ef = genEfilt(tau,fBins);%exponential filter
ef_init = ef;

samples_a  = cell(1,nsweeps);
samples_s = cell(1,nsweeps);
samples_pr = cell(1,nsweeps);
samples_tau = cell(1,nsweeps);
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

efs = cell(1,length(Tguess));

pr = b*ones(1,nBins); %initial calcium is set to baseline 

N = length(sti); %number of spikes in spiketrain

%initial logC - compute likelihood initially completely - updates to likelihood will be local
logC = sum(-(pr-trace).^2); 

m = p_spike*nBins;

logC_ = logC;
for i = 1:length(Tguess)
    efs{i} = ef;
    tmpi = Tguess(i); 
    sti_ = [sti tmpi];
    %must add spike to each trial (at mean location or sampled -- more appropriate if sampled)
    pr_ = pr;
    ati_ = ati;
    a_init = max(trace(tmpi)/A,a_min);
    [sti_, pr_, logC_] = addSpike(sti,pr,logC_,efs{i},a_init,tau,trace,tmpi, N+1, Dt, A); %adds all trials' spikes at same time
    taus{i} = tau;
    ati_ = [ati_ a_init];
    ati = ati_;
    sti = sti_;
    pr = pr_;
    N = length(sti); %number of spikes in spiketrain
end
logC = logC_;

sti_= sti;
logC_= logC;
N=length(sti);

%% loop over sweeps to generate samples
addMoves = [0 0]; %first elem is number successful, second is number total
dropMoves = [0 0];
timeMoves = [0 0];
ampMoves = [0 0];
tauMoves = [0 0];
for i = 1:nsweeps
%     
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
            [si_, pr_, logC_] = removeSpike(si,pr,logC,efs{ni},ai(ni),taus{ni},trace,tmpi,ni, Dt, A);
            [si_, pr_, logC_] = addSpike(si_,pr_,logC_,efs{ni},ai(ni),taus{ni},trace,tmpi_,ni, Dt, A);

            %accept or reject
            %for prior: (1) use ratio or(2) set prior to 1.
            prior_ratio = 1;

            ratio = exp(sum((1/(2*NoiseVar))*(logC_-logC)))*prior_ratio;
            if ratio>1 %accept
                si = si_;
                pr = pr_;
                logC = logC_;
                timeMoves = timeMoves + [1 1];
                proposalVar = proposalVar + 2*.1*rand*proposalVar/sqrt(i);
            elseif rand<ratio %accept
                si = si_;
                pr = pr_;
                logC = logC_;
                timeMoves = timeMoves + [1 1];
                proposalVar = proposalVar + 2*.1*rand*proposalVar/sqrt(i);
            else
                %reject - do nothing
                proposalVar = proposalVar - .1*rand*proposalVar/sqrt(i);
                timeMoves = timeMoves + [0 1];
            end
        end
        sti = si;
    end

    
    % update amplitude of each burst
    for ii = 1:10
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
            [si_, pr_, logC_] = removeSpike(si,pr,logC,efs{ni},ai(ni),taus{ni},trace,si(ni),ni, Dt, A);
            [si_, pr_, logC_] = addSpike(si_,pr_,logC_,efs{ni},tmp_a_,taus{ni},trace,si(ni),ni, Dt, A);

            ai_ = ai;
            ai_(ni) = tmp_a_;

            %accept or reject - include a prior?
            prior_ratio = 1;
            ratio = exp(sum((1/(2*NoiseVar))*(logC_-logC)))*prior_ratio;
            if ratio>1 %accept
                ai = ai_;
                si = si_;
                pr = pr_;
                logC = logC_;
                ampMoves = ampMoves + [1 1];
                a_std = a_std + 2*.1*rand*a_std/sqrt(i);
            elseif rand<ratio %accept
                ai = ai_;
                si = si_;
                pr = pr_;
                logC = logC_;
                ampMoves = ampMoves + [1 1];
                a_std = a_std + 2*.1*rand*a_std/sqrt(i);
            else
                %reject - do nothing
                a_std = a_std - .1*rand*a_std/sqrt(i);
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
        [pr_, logC_] = remove_base(pr,logC,tmp_b,trace,A);   
        [pr_, logC_] = add_base(pr_,logC_,tmp_b_,trace,A);

        %accept or reject - include a prior?
        prior_ratio = 1;
        ratio = exp(sum((1/(2*NoiseVar))*(logC_-logC)))*prior_ratio;
        if ratio>1 %accept
            baseline = tmp_b_;
            pr = pr_;
            logC = logC_;
            b_std = b_std + 2*.1*rand*b_std/sqrt(i);
        elseif rand<ratio %accept
            baseline = tmp_b_;
            pr = pr_;
            logC = logC_;
            b_std = b_std + 2*.1*rand*b_std/sqrt(i);
        else
            b_std = b_std - .1*rand*b_std/sqrt(i);
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
%             [~,tmpi] = max((pr-trace).^2); % this problem isn't a legit
%             proposal
            %dont add if we have too many bursts or the proposed new location
            %is too close to another one
            if ~(any(abs(tmpi-sti)<exclusion_bound) || N >= maxNbursts)
                sti_ = [sti tmpi];
                %must add burst to each trial (at mean location or sampled -- more appropriate if sampled, but make sure no trial's burst violates exclusion)
                logC_ = logC;
                pr_ = pr;
                ati_ = ati;
                a_init = max(trace(max(1,floor(tmpi)))/A - baseline + a_std*randn,a_min);%propose an initial amplitude for it
                [si_, pr_, logC_] = addSpike(sti,pr,logC_,ef_init,a_init,tau,trace,tmpi, N+1, Dt, A); %adds all trials' bursts at same time
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
                ratio = exp(sum((1./(2*NoiseVar)).*(logC_-logC)))*(rprob/fprob)*(m(1)/(nBins(1)-m(1))); %posterior times reverse prob/forward prob
                if (ratio>1)||(ratio>rand) %accept
                    ati = ati_;
                    sti = sti_;
                    pr = pr_;
                    taus{N+1} = tau;
                    efs{N+1} = ef_init;
                    logC = logC_;
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
                logC_ = logC;
                pr_ = pr;
                ati_ = ati;
                %always remove the ith burst (the ith burst of each trial is linked)                     
                [si_, pr_, logC_] = removeSpike(sti,pr,logC_,efs{tmpi},ati(tmpi),taus{tmpi},trace,sti(tmpi),tmpi, Dt, A);
                sti_ = si_;
                ati_(tmpi) = [];

                %reverse probability
                rprob = 1/nBins(1);

                %compute forward prob
                fprob = 1/N;

                %accept or reject
                %posterior times reverse prob/forward prob
                ratio = exp(sum((1./(2*NoiseVar)).*(logC_-logC)))*(rprob/fprob)*((nBins(1)-m(1))/m(1)); 
                if (ratio>1)||(ratio>rand)%accept
                    ati = ati_;
                    sti = sti_;
                    pr = pr_;
                    taus(tmpi) = [];
                    efs(tmpi) = [];
                    logC = logC_;
                    dropMoves = dropMoves + [1 1]; 
                else
                    %reject - do nothing
                    dropMoves = dropMoves + [0 1];
                end
                N = length(sti);
            end
        end
    end
    
    
    %% this is the section that updates tau
    % update tau (via random walk sampling)
    for ii = 1:1
        for ni = 1:N 
            % update both tau values
            tau_ = taus{ni};
            tau_(1) = tau_(1)+(tau1_std*randn); %with bouncing off min and max
%             count1 = 0;
%             count2 = 0;
           while tau_(1)>tau(2) || tau_(1)<tau_min
                if tau_(1)<tau_min
                    tau_(1) = tau_min+(tau_min-tau_(1));
                elseif tau_(1)>tau(2)
                    tau_(1) = tau(2)-(tau_(1)-tau(2));
                end
            end 

            ef_ = genEfilt(tau_,fBins);%exponential filter

            %remove all old bumps and replace them with new bumps    
            logC_ = logC;
            pr_ = pr;
            [~, pr_, logC_] = removeSpike(sti,pr_,logC_,efs{ni},ati(ni),taus{ni},trace,sti(ni),ni, Dt, A);
            [~, pr_, logC_] = addSpike(sti,pr_,logC_,ef_,ati(ni),tau_,trace,sti(ni),ni, Dt, A);

            %accept or reject
            prior_ratio = 1;
%             prior_ratio = (gampdf(tau_(1),1.5,1))/(gampdf(tau(1),1.5,1));
            ratio = exp(sum(sum((1./(2*NoiseVar)).*(logC_-logC))))*prior_ratio;
            if ratio>1 %accept
                pr = pr_;
                logC = logC_;
                taus{ni} = tau_;
                efs{ni} = ef_;
                tauMoves = tauMoves + [1 1];
%                 tau1_std = tau1_std + 2*.1*rand*tau1_std/sqrt(i);
            elseif rand<ratio %accept
                pr = pr_;
                logC = logC_;
                taus{ni} = tau_;
                efs{ni} = ef_;
                tauMoves = tauMoves + [1 1];
%                 tau1_std = tau1_std + 2*.1*rand*tau1_std/sqrt(i);
            else
                %reject - do nothing
%                 tau1_std = tau1_std - .1*rand*tau1_std/sqrt(i);
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
%             count1 = 0;
%             count2 = 0;
            while tau_(2)>tau_max || tau_(2)<tau_(1)
                if tau_(2)<tau_(1)
                    tau_(2) = tau_(1)+(tau_(1)-tau_(2));
                elseif tau_(2)>tau_max
                    tau_(2) = tau_max-(tau_(2)-tau_max);
                end
            end   

            ef_ = genEfilt(tau_,fBins);%exponential filter

            %remove all old bumps and replace them with new bumps    
            logC_ = logC;
            pr_ = pr;
            [~, pr_, logC_] = removeSpike(sti,pr_,logC_,efs{ni},ati(ni),taus{ni},trace,sti(ni),ni, Dt, A);
            [~, pr_, logC_] = addSpike(sti,pr_,logC_,ef_,ati(ni),tau_,trace,sti(ni),ni, Dt, A);

            %accept or reject
            prior_ratio = 1;
%             prior_ratio = gampdf(tau_(2),12,1)/gampdf(tau(2),12,1);
            ratio = exp(sum(sum((1./(2*NoiseVar)).*(logC_-logC))))*prior_ratio;
            if ratio>1 %accept
                pr = pr_;
                logC = logC_;
                taus{ni} = tau_;
                efs{ni} = ef_;
                tauMoves = tauMoves + [1 1];
%                 tau2_std = tau2_std + 2*.1*rand*tau2_std/sqrt(i);
            elseif rand<ratio %accept
                pr = pr_;
                logC = logC_;
                taus{ni} = tau_;
                efs{ni} = ef_;
                tauMoves = tauMoves + [1 1];
%                 tau2_std = tau2_std + 2*.1*rand*tau2_std/sqrt(i);
            else
                %reject - do nothing
%                 tau2_std = tau2_std - .1*rand*tau2_std/sqrt(i);
                tauMoves = tauMoves + [0 1];
            end
        end
    end
    
    
    
    % re-estimate the noise variance
    if ~isempty(sti)
        df = (numel(pr)); %DOF (possibly numel(ci(ti,:))-1)
        NoiseVar = sum((pr-trace).^2)/df; %ML - DOF, init, baseline and each burst amplitude

        A_samp = 0.5 * (df);
        B_samp = 1/(0.5 * df * NoiseVar);
        NoiseVar = 1/gamrnd(A_samp,B_samp); %this could be inf but it shouldn't be
    end
    

    %store things
    N_sto = [N_sto N];
    samples_a{i} = ati; %trial amplitudes
    samples_b{i} = baseline; %trial baselines
    samples_s{i} = sti; %shared bursts
    samples_pr{i} = pr; %save calcium traces
    samples_tau{i} = taus; %save tau values
    %store overall logliklihood as well
%     if abs(sum(logC)-sum(sum(-(pr)-cell2mat(trace)).^2)))>1
%         figure(90)
%         subplot(121)
%         plot(cell2mat(samples_c{i-1})')``
%         subplot(122)
%         plot(cell2mat(pr)')
%         keyboard
%     end

    objective = [objective logC];
%     figure(10);
%     plot(ci{1});hold on;
%     plot(CaF{1},'r');hold off
%     if sum(ismember(indreporti,i))
%         fprintf([num2str(indreport(ismember(indreporti,i)),2),', '])
%     end
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
mcmc.noisevar = NoiseVar;

trials.amp=samples_a;
trials.base=samples_b;
% trials.curves=samples_pr;
trials.tau=samples_tau;
trials.obj = -objective;
trials.times = samples_s;

params.NoiseVar = NoiseVar_init; 
params.proposalVar = proposalVar;
params.nsweeps = nsweeps;
params.tau1_std = tau1_std; %proposal variance of tau parameters
params.tau2_std = tau2_std; %proposal variance of tau parameters
params.tau_min = tau_min;
params.tau_max = tau_max;

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


params.adddrop = adddrop;
params.maxNbursts = maxNbursts;

% 
% disp('Below are the moves that were done')
% display(['time: ' num2str(timeMoves(1)/timeMoves(2))])
% display(['add: ' num2str(addMoves(1)/addMoves(2))])
% display(['drop: ' num2str(dropMoves(1)/dropMoves(2))])
% display(['amplitude: ' num2str(ampMoves(1)/ampMoves(2))])
% display(['tau: ' num2str(tauMoves(1)/tauMoves(2))])




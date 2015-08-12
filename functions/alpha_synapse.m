function current = alpha_synapse(t, onset, tau, gmax)
%ALPHA_SYNAPSE
%
% based on http://kaplan.bio.upenn.edu/nia2006/nia/NEUROLAB/help/ALPHASYN.HTM
%
% arguments:
%
%      t: range of times to consider (ms). for example, 0:.01:2
%  onset: time (ms) that change in postsynaptic conductance begins; e.g. 0.5
%    tau: time (ms) to peak of conductance change; e.g. 0.1 
%   gmax: peak conductance change, in micro-Siemens (micro-mho); e.g. 2
%      e: reversal potential for synaptic current, in mV; e.g. -15
%      v: resting potential (in mV); e.g. 0
% 
% returns:
%
% synaptic current, in nA, over the time domain t
%

onset = onset(:);
n_onsets = length(onset);

t_matrix = t(ones(n_onsets,1),:);

onset_matrix = onset(:,ones(length(t),1));
valid_t_ix = t_matrix > onset_matrix;
t_gt_onset = t_matrix(valid_t_ix);
onset_matrix = onset_matrix(valid_t_ix);

%t_gt_onset = t(t > onset);

current = zeros(size(t_matrix));
current(valid_t_ix) = gmax .* (t_gt_onset - onset_matrix)./tau .* exp(-(t_gt_onset - onset_matrix - tau)./tau);

% current = g .* (v - e);


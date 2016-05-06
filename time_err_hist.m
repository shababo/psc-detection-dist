function errs = time_err_hist(dirname,match_string,param_variables)

dirinfo = dir(dirname);

a_mins = param_variables.a_mins;
p_spikes = param_variables.p_spikes;
num_samps = param_variables.num_samps;
roc = zeros(length(a_mins),length(p_spikes),length(num_samps),2);

errs = [];
for i = 1:length(dirinfo)
    
    dirinfo(i).name

    if regexpi(dirinfo(i).name,match_string)
        
        disp('matched')
        
        load([dirname '/' dirinfo(i).name],'params')
        
        a_min_i = find(a_mins == params.a_min);
        p_spike_i = find(p_spikes == params.p_spike);
        num_samps_i = find(num_samps == params.num_sweeps);
        
        if isempty(a_min_i) || isempty(p_spike_i) || isempty(num_samps_i) || num_samps(num_samps_i) < 500
            continue
        end
        
        load([dirname '/' dirinfo(i).name])
        if ~exist('timing_score','var') || isstr(results)
            disp(['not scored: ' dirinfo(i).name])
            continue
        end
        
     
        for j = 1:length(results)
            
            map_ind = results(j).map_ind;
            num_est_events = length(results(j).trials.times{map_ind});
            
            
%             if length(num_est_events) - length(timing_score(j).correct_inds) > 40
%                 disp('too many fp')
%                 continue
%             else
%                 num_est_events - length(timing_score(j).correct_inds)
%             end
            errs = [errs timing_score(j).correct_err];
            
        end
    end
end

figure;
histogram(errs)
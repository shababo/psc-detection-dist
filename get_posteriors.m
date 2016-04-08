function posteriors = get_posteriors(filename)

load(filename)

posteriors = struct();

for i = 1:length(results)
    
    % get amplitudes
    posteriors(i).amps = [];
    posteriors(i).num_events = [];
    posteriors(i).tau1 = [];
    posteriors(i).tau2 = [];
    posteriors(i).times = [];
    
    for j = 1:length(results(i).trials.amp)
        
         posteriors(i).amps = [posteriors(i).amps results(i).trials.amp{j}];
         posteriors(i).num_events = [posteriors(i).num_events length(results(i).trials.amp{j})];
         
         for k = 1:length(results(i).trials.amp{j})
             posteriors(i).tau1 = [posteriors(i).tau1 results(i).trials.tau{j}{k}(1)];
             posteriors(i).tau2 = [posteriors(i).tau2 results(i).trials.tau{j}{k}(2)];
         end
         
         posteriors(i).times = [posteriors(i).times results(i).trials.times{j}];
         
    end
    
end
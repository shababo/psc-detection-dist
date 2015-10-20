function errs = times_err_hist_deconv(filename)

load(filename)

errs = [];

for thresh_i = 1:length(threshold)
    for min_spacing_i = 1:length(min_spacing)
       
        if ~exist('timing_score','var')
            disp(['not scored:'])
            continue
        end
        
        for j = 1:length(event_times{thresh_i,min_spacing_i})
            
%             if length(event_times{thresh_i,min_spacing_i}{j}) - length(timing_score(thresh_i,min_spacing_i,j).correct_err) > 40
%                 disp('too many false pos')
%                 continue
%             else
%                 disp('not too many')
%                 length(event_times{thresh_i,min_spacing_i}{j}) - length(timing_score(thresh_i,min_spacing_i,j).correct_err)
%             end
            errs = [errs timing_score(thresh_i,min_spacing_i,j).correct_err];
            
            
        end
    end
end

figure
histogram(errs,'normalization','pdf')
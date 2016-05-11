function roc = compute_roc_deconv(filename)

load(filename)

roc = zeros(size(event_times,1),size(event_times,2),2);

for thresh_i = 1:size(event_times,1)
    for min_spacing_i = 1:size(event_times,2)
       
        if ~exist('timing_score','var')
            disp(['not scored: '])
            continue
        end
        
        for j = 1:length(event_times{thresh_i,min_spacing_i})
            
            num_est_events = length(event_times{thresh_i,min_spacing_i}{j});
            roc(thresh_i,min_spacing_i,1) = roc(thresh_i,min_spacing_i,1) + num_est_events - length(timing_score(thresh_i,min_spacing_i,j).correct_inds);
            roc(thresh_i,min_spacing_i,2) = roc(thresh_i,min_spacing_i,2) + length(timing_score(thresh_i,min_spacing_i,j).correct_inds);
            
        end
    end
end
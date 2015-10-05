function [correct_est_inds, correct_est_errs] = score_detection_trace(est_times, true_times, tolerance)

correct_est_inds = [];
correct_est_errs = [];

est_times_full = est_times;

for i = 1:length(true_times)
    
    [min_offset, min_offset_ind] = min(abs(est_times - true_times(i)));
    if min_offset < tolerance
        correct_est_inds = [correct_est_inds find(est_times_full == est_times(min_offset_ind))];
        correct_est_errs = [correct_est_errs min_offset];
        est_times(min_offset_ind) = [];
    end
end


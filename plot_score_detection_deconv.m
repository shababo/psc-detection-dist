function plot_score_detection_deconv(est_event_times,true_event_times,timing_score)


for i = 1:length(est_event_times)
    i
    est_times = est_event_times{i};
    scatter(true_event_times{i},-i*ones(size(true_event_times{i})),'b*')
    hold on
    correct_times = est_times(timing_score(i).correct_inds);
    false_pos_times = est_times(setdiff(1:length(est_times),timing_score(i).correct_inds));
    scatter(correct_times,-i*ones(size(correct_times)),'go')
    hold on
    scatter(false_pos_times,-i*ones(size(false_pos_times)),'rx')
end
hold off
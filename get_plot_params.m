function plot_params = get_plot_params

plot_params.min_tau1 = 0;
plot_params.min_tau2 = 0;
plot_params.max_tau1 = Inf;
plot_params.max_tau2 = Inf;

plot_params.min_amp = 5;
plot_params.min_time = 1;
plot_params.max_time = Inf;%.015*20000;
plot_params.hot_time_min = .0473;
plot_params.hot_time_max = .0494;
% plot_params.hot_time_min = .1519;
% plot_params.hot_time_max = .1526;
% plot_params.hot_time_min = .052;
% plot_params.hot_time_max = .05275;
% plot_params.hot_time_min = .2521;
% plot_params.hot_time_max = .2531;
% plot_params.hot_time_min = .0047;
% plot_params.hot_time_max = .007;
plot_params.hot_amp_min = 0;
% plot_params.hot_amp_min = 20;
plot_params.hot_amp_max = Inf;
plot_params.hot_min_tau2 = 0;%50;


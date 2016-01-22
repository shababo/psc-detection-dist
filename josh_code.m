mean_AR2 = zeros(length(noise_vals),length(simTrace));
mean_AR0 = zeros(length(noise_vals),length(simTrace));
fs = 2e4;
for ni = 1
    load('data/random_snr_3_rep_1.mat') 
    figure(9)
    subplot(611)
    plot((1:N)/fs,-mc_sim,'k')
    title('simulated event trace')
    ylim([-15 0])
    xlim([2*.2e4 2*.5e4]./fs)
    set(gca,'XTick',[])
    set(gca,'xColor','w')
    box off
    subplot(612)
    plot((1:N)/fs,-er(1:N),'k')
    title('error process')
    ylim([-7.5 7.5])
    xlim([2*.2e4 2*.5e4]./fs)
    set(gca,'XTick',[])
    set(gca,'xColor','w')
    box off
    subplot(613)
    plot((1:N)/fs,-simTrace,'k')
    title('observed trace (events + noise)')
    ylim([-15 0])
    xlim([2*.2e4 2*.5e4]./fs)
    set(gca,'XTick',[])
    set(gca,'xColor','w')
    box off
    subplot(614)
    mean_AR2(ni,:) = PlotMCMC(results1,simTrace);
    hold on; plot(trace_info_sim.trials.times{1}/fs,-1*ones(n_event,1),'.b','MarkerFaceColor','b','markersize',5)
    title('inference with AR(2) noise model')
    ylim([-15 0])
    xlim([2*.2e4 2*.5e4]./fs)
    set(gca,'XTick',[])
    set(gca,'xColor','w')
    box off
    subplot(615)
    mean_AR0(ni,:) = PlotMCMC(results2,simTrace);
    hold on; plot(trace_info_sim.trials.times{1}/fs,-1*ones(n_event,1),'.b','MarkerFaceColor','b','markersize',5)
    title('inference with AR(0) noise model')
    ylim([-15 0])
    xlim([2*.2e4 2*.5e4]./fs)
    set(gca,'XTick',[])
    set(gca,'xColor','w')
    box off
    subplot(616)
    plot(trace_info_sim.trials.times{1}/fs,-trace_info_sim.trials.amp{1},'.k')
    hold on
    for i = 400:length(results1.trials.amp)
        plot(results1.trials.times{i}/fs,-results1.trials.amp{i},'.r')
    end
    plot(trace_info_sim.trials.times{1}/fs,-trace_info_sim.trials.amp{1},'.k')
    hold off
    xlabel('time (s)')
    ylabel('amp')
    title('posterior joint samples of event time and amplitude AR(2) model')
    ylim([-15 0])
    xlim([2*.2e4 2*.5e4]./fs)
    set(gcf,'position',[185   289   731   470]);
    box off
end
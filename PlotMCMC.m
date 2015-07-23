
%% multiple samples for same trace
figure(1)
burnIn = round(2/5*length(results.trials.tau));
nBins = length(trace);

plot(trace,'ko'); 
modelledTrace = [];
hold on
for i = burnIn:10:length(results.trials.tau)
    modelledTrace = [modelledTrace; plot_curve(results, trace, i)];
    plot_curve(results, trace, i, 'r');
%     plot(trials.curves{i},'r');
    xlabel('time (ms)')
    axis tight
end
MTm = mean(modelledTrace);
MTs = std(modelledTrace);%estimation noise
plot(MTm,'g');
hold off
axis tight
hold off
axis tight

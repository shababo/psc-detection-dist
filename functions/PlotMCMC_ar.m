
%% multiple samples for same trace
figure(1)
burnIn = round(2/5*length(trials.tau));
nBins = length(trace);
%ML - DOF, init, baseline and each burst amplitude
SSE=sum((trace-trials.curves{end}).^2);
n=numel(trials.curves{end})-(2+mcmc.N_sto(end)); 
final_var = SSE/n;%standard error
%
final_std = sqrt(final_var);
plot(trace,'ko'); 
modelledTrace = [];
hold on
for i = burnIn:10:length(trials.curves)
    modelledTrace = [modelledTrace;trials.curves{i}];
    plot(trials.curves{i},'r');
    xlabel('time (ms)')
    axis tight
end
MTm = mean(modelledTrace);
MTs = std(modelledTrace);%estimation noise
plot(MTm,'g');
%Liams idea: posterior uncertaubnty which is the sum of uncertainty in
%Ca signal + inherent noise of calcium signal
plot(MTm+2*(final_std+MTs),'r--'); 
plot(MTm-2*(final_std+MTs),'r--'); 
hold off
axis tight
hold off
axis tight

function [newSpikeTrain, newCalcium, newLL, tmp] = addSpike(oldSpikeTrain,oldCalcium,oldLL,filter,amp,tau,obsCalcium,timeToAdd,indx,Dt,A)

    tau_h = tau(1);
    tau_d = tau(2);

    ef_h = filter{1};
    ef_d = filter{2};

    %     newSpikeTrain = [oldSpikeTrain timeToAdd]; %possibly inefficient, change if problematic (only likely to be a problem for large numbers of spikes)
    newSpikeTrain = [oldSpikeTrain(1:indx-1) timeToAdd oldSpikeTrain(indx:end)]; %adds in place rather than at end

    %use infinite precision to scale the precomputed FIR approximation to the calcium transient
    wk_h = amp*A*exp((timeToAdd - Dt*ceil(timeToAdd/Dt))/tau_h);
    wk_d = amp*A*exp((timeToAdd - Dt*ceil(timeToAdd/Dt))/tau_d);    
        
    %%%%%%%%%%%%%%%%%
    %handle ef_h first
    newCalcium = oldCalcium;
    tmp = 1 + (floor(timeToAdd):min((length(ef_h)+floor(timeToAdd)-1),length(newCalcium)-1));
    newCalcium(tmp) = newCalcium(tmp) + wk_h*ef_h(1:length(tmp));

    %if you really want to, ef*ef' could be precomputed and passed in
%     relevantResidual = obsCalcium(tmp)-oldCalcium(tmp);
%     newLL = oldLL - ( wk_h^2*norm(ef_h(1:length(tmp)))^2 - 2*relevantResidual*(wk_h*ef_h(1:length(tmp))'));
    newLL = oldLL;
    newLL(tmp) = obsCalcium(tmp) - newCalcium(tmp);
% % %     display(['add1 indx:' num2str(timeToAdd) ' ' num2str(newLL-oldLL)])
% % %     keyboard
    oldCalcium = newCalcium;
%     oldLL = newLL;
    %%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%
    %handle ef_d next
    tmp = 1 + (floor(timeToAdd):min((length(ef_d)+floor(timeToAdd)-1),length(newCalcium)-1));
    newCalcium(tmp) = newCalcium(tmp) + wk_d*ef_d(1:length(tmp));

    %if you really want to, ef*ef' could be precomputed and passed in
%     relevantResidual = obsCalcium(tmp)-oldCalcium(tmp);
%     newLL = oldLL - ( wk_d^2*norm(ef_d(1:length(tmp)))^2 - 2*relevantResidual*(wk_d*ef_d(1:length(tmp))'));
% % %     display(['add2 indx:' num2str(timeToAdd) ' ' num2str(newLL-oldLL)])
% % %     keyboard
    newLL(tmp) = obsCalcium(tmp) - newCalcium(tmp);
    
    %%%%%%%%%%%%%%%%%
    %separately handle LL
%     newLL = -sum((newCalcium-obsCalcium).^2);





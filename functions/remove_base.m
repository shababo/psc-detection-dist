function [newCalcium, newLL] = remove_base(oldCalcium,oldLL,amp,obsCalcium,A)
    
    newCalcium = oldCalcium;
    
    newCalcium = newCalcium - amp*A*ones(1,length(oldCalcium));

%     %% must recompute error
%     oldResidual = obsCalcium - oldCalcium;
%     newResidual = obsCalcium - newCalcium;
%     newLL = oldLL - sum(-(oldResidual).^2) + sum(-(newResidual).^2);
   
    newLL = -sum((newCalcium-obsCalcium).^2);
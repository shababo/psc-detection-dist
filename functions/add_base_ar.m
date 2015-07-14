function [newCalcium, newLL] = add_base(oldCalcium,oldLL,amp,obsCalcium,A)
    
newCalcium = oldCalcium;

newCalcium = newCalcium + amp*A*ones(1,length(oldCalcium));
% % must recompute error
% oldResidual = obsCalcium - oldCalcium;
% newResidual = obsCalcium - newCalcium;
% % oldLL
% % sum(-(oldResidual).^2)
% % sum(-(newResidual).^2)
% newLL = oldLL - sum(-(oldResidual).^2) + sum(-(newResidual).^2);

newLL = obsCalcium-newCalcium;

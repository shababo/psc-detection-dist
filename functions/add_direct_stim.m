function [new_trace, diffY] = add_direct_stim(trace,stim_amp,diffY,shape)

new_trace = trace + stim_amp*shape;

diffY = diffY - stim_amp*shape;
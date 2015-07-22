function [new_trace, logC] = add_line_noise(trace,phi,a_s,theta,obs)

new_trace = trace + a_s*sin(((1:length(trace))*theta + phi)*2*pi);

logC = -sum((new_trace - obs).^2);
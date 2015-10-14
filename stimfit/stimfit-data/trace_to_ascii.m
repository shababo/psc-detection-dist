function trace_to_ascii(filename,trace,tag,dt)

dlmwrite(filename,tag,'')
dlmwrite(filename,length(trace),'-append','delimiter','')
trace = reshape(trace,[length(trace) 1]);
time = ((0:(length(trace) - 1))*dt)';
dlmwrite(filename,[time trace],'-append','delimiter',' ');
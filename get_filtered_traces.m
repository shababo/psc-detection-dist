function filtered_traces = get_filtered_traces(results)

filtered_traces = zeros(length(results),length(results(1).filtered_trace));

for i = 1:length(results)
    filtered_traces(i,:) = results(i).filtered_trace;
end
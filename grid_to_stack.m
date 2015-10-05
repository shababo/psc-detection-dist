function traces = grid_to_stack(trace_array)

num_rows = size(trace_array,1);
num_cols = size(trace_array,2);

num_traces = 0;
for i = 1:num_rows
    for j = 1:num_cols
        num_traces = num_traces + size(trace_array{i,j},2);
    end
end

duration = size(trace_array{1,1},1);

traces = zeros(num_traces,duration);

count = 1;
for i = 1:num_rows
    for j = 1:num_cols
        offset = size(trace_array{i,j},2) - 1;
        traces(count:count+offset,:) = trace_array{i,j}';
        count = count + offset + 1;
    end
end
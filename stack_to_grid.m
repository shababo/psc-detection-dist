function trace_array_new = stack_to_grid(traces,trace_array)

num_rows = size(trace_array,1);
num_cols = size(trace_array,2);

trace_array_new = cell(num_rows,num_cols);

ind = 1;
for i = 1:num_rows
    for j = 1:num_cols
        
        offset = size(trace_array{i,j},2)-1;
        trace_array_new{i,j} = traces(ind:ind+offset,:)';
        ind = ind+offset+1;
        
    end
end
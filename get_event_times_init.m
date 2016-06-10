function results_out = get_event_times_init(results,trace_length,do_bit_vec,min_size)

if do_bit_vec
    results_out = zeros(length(results),trace_length);
else
    results_out = cell(length(results),1);
end

for i = 1:length(results)
    if do_bit_vec
        for j = 1:length(results(i).event_times_init)
            if results(i).event_sizes_init >= min_size
                results_out(i,results(i).event_times_init(j)) = 1;
            end

        end
    else
        results_out{i} = [];
        for j = 1:length(results(i).event_times_init)
            if results(i).event_sizes_init >= min_size
                results_out(i,results(i).event_times_init(j)) = 1;
                results_out{i} = [results_out{i} results(i).event_times_init(j)];
            end

        end
        
    end
end
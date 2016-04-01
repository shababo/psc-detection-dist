function results_out = get_event_times_init(results,trace_length,do_bit_vec)

if do_bit_vec
    results_out = zeros(length(results),trace_length);
else
    results_out = cell(length(results),1);
end

for i = 1:length(results)
    if do_bit_vec
        for j = 1:length(results(i).event_times_init)
            results_out(i,results(i).event_times_init) = 1;
        end
    else
        results_out{i} = results(i).event_times_init;
    end
end
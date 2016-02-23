function mean_n_grid = get_mean_N_grid(results_grid,burn_in)

mean_n_grid = cell(size(results_grid));

for i = 1:size(results_grid,1)
    for j = 1:size(results_grid,2)
        mean_n_grid{i,j} = zeros(length(results_grid{i,j}),1);
        for k = 1:length(results_grid{i,j})
            mean_n_grid{i,j}(k) = mean(results_grid{i,j}(k).mcmc.N_sto(burn_in:end));
        end
    end
end
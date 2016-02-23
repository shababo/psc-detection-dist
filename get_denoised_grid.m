function denoised_grid = get_denoised_grid(results_grid,params,T)

denoised_grid = cell(size(results_grid));

for i = 1:size(results_grid,1)
    for j = 1:size(results_grid,2)
        denoised_grid{i,j} = zeros(length(results_grid{i,j}),T);
        for k = 1:length(results_grid{i,j})
            denoised_grid{i,j}(k,:) = get_map_trace(results_grid{i,j}(k),params,T);
        end
    end
end
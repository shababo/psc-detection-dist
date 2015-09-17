function plot_obj(results_file,linespec)

load(results_file)

colors = lines(length(results));

for i = 1:length(results)
    
    plot(results(i).trials.obj,linespec,'Color',colors(i,:))
    hold on
end

hold off
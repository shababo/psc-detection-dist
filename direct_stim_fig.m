load('data/direct-stim-work-2P.mat')
%%
figure
subplot(121); 
plot(template_traces([1:2 4:10],:)')
hold on;

bar_corner_time = 1500;
bar_corner_y = -150;
bar_limits = [300 25];
        plot([bar_corner_time; bar_corner_time], bar_corner_y + [0; bar_limits(2)], '-k',  bar_corner_time + [0; bar_limits(1)], [bar_corner_y; bar_corner_y], '-k', 'LineWidth', 2)
        text(bar_corner_time - bar_limits(1)/2,bar_corner_y + bar_limits(2)/2, [num2str(bar_limits(2)) ' pA'], 'HorizontalAlignment','right')
        text(bar_corner_time + bar_limits(1)/2,bar_corner_y - bar_limits(2)/2, [num2str(bar_limits(1)/20) ' ms'], 'HorizontalAlignment','center')
axis tight
axis off
title('Direct Optical Current')

subplot(122)
plot(norm_template_traces([1:2 4:10],:)')
title('Normalized Optical Currents')
axis tight
axis off
tau1 = 1:2:20;
tau2 = 10:20:200;


figure
c = 1;
rise_times = zeros(length(tau1),length(tau2));
decay_times = zeros(length(tau1),length(tau2));
max_times = zeros(length(tau1),length(tau2));
max_val = zeros(length(tau1),length(tau2));
for i = 1:length(tau1)
    for j = 1:length(tau2)
        subplot(length(tau1),length(tau2),c)
        this_curve = zeros(1,800);
    ef = genEfilt_ar([tau1(i) tau2(j)],400);
    [~, this_curve, ~] = addSpike_ar(1000,...
                                            this_curve, 0, ef, ...
                                            1,...
                                            [tau1(i) tau2(i)],...
                                            zeros(size(this_curve)), ...
                                            1, ...
                                            2, 1, 1);
                                        c = c+1;
                                        
                                        [max_val(i,j),max_times(i,j)] = max(this_curve);
                                        rise_times(i,j) = find(this_curve > .67*max(this_curve),1,'first');
                                        
                                        decay_times(i,j) = find(this_curve(rise_times(i,j):end) < .33*max(this_curve),1,'first') + rise_times(i,j) - max_times(i,j);
    
    
    
        plot(this_curve,'LineWidth',2)
        hold on
        scatter(rise_times(i,j),.67*max(this_curve),'g')
        hold on
        scatter(decay_times(i,j),.33*max(this_curve),'r')
        hold off
                                        axis off
                                        
                                     
                                        
                                        
    end
end

   figure;
                                        subplot(1,4,1)
                                        imagesc(rise_times)
                                        subplot(1,4,2)
                                        imagesc(decay_times)
                                        subplot(1,4,3)
                                        imagesc(max_times)
                                        subplot(1,4,4)
                                        imagesc(max_val)
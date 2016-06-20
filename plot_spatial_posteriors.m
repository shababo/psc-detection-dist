function plot_spatial_posteriors(spatial_posteriors)

% spatial_posteriors = struct();


figure(9991)
subplot1(size(spatial_posteriors,1), size(spatial_posteriors,2), 'Gap', [.005 .005], 'XTickL', 'Margin', 'YTickL', 'Margin');
figure(9992)
subplot1(size(spatial_posteriors,1), size(spatial_posteriors,2), 'Gap', [.005 .005], 'XTickL', 'Margin', 'YTickL', 'Margin');
% figure(9993)
% subplot1(size(spatial_posteriors,1), size(spatial_posteriors,2), 'Gap', [.005 .005], 'XTickL', 'Margin', 'YTickL', 'Margin');
figure(99944)
subplot1(size(spatial_posteriors,1), size(spatial_posteriors,2), 'Gap', [.005 .005], 'XTickL', 'Margin', 'YTickL', 'Margin');
figure(9995)
subplot1(size(spatial_posteriors,1), size(spatial_posteriors,2), 'Gap', [.005 .005], 'XTickL', 'Margin', 'YTickL', 'Margin');

for i = 1:size(spatial_posteriors,1)
    for j = 1:size(spatial_posteriors,2)
        
%         if ~isempty(spatial_posteriors{i,j})
%             
%             spatial_posteriors(i,j).amps = [];
%             spatial_posteriors(i,j).num_events = [];
%             spatial_posteriors(i,j).tau1 = [];
%             spatial_posteriors(i,j).tau2 = [];
%             spatial_posteriors(i,j).times = [];
%             
%             for ii = -neighborhood_size:neighborhood_size
%                 for jj = -neighborhood_size:neighborhood_size
%                      if ~(i+ii < 1 || i+ii > size(spatial_posteriors,1) || j+jj < 1 || j+jj > size(spatial_posteriors,2))
%                          
%                          133
%                          spatial_posteriors(i,j).amps = [spatial_posteriors(i,j).amps [spatial_posteriors{i+ii,j+jj}.amps]];
%                          spatial_posteriors(i,j).num_events = [spatial_posteriors(i,j).num_events [spatial_posteriors{i+ii,j+jj}.num_events]];
%                          spatial_posteriors(i,j).tau1 = [spatial_posteriors(i,j).tau1 [spatial_posteriors{i+ii,j+jj}.tau1]];
%                          spatial_posteriors(i,j).tau2 = [spatial_posteriors(i,j).tau2 [spatial_posteriors{i+ii,j+jj}.tau2]];
%                          spatial_posteriors(i,j).times = [spatial_posteriors(i,j).times [spatial_posteriors{i+ii,j+jj}.times]];
%                          
%                      end
%                 end
%             end

            
            (i-1)*size(spatial_posteriors,2) + j

                        
            figure(9991)
            subplot1((i-1)*size(spatial_posteriors,2) + j);
            histogram(spatial_posteriors(i,j).amps,[20:3:300],'Normalization','countdensity')
            axis off
            ylim([0 40]*10*3)
            xlim([20 300])
            
            figure(9992)
            subplot1((i-1)*size(spatial_posteriors,2) + j);
            histogram(spatial_posteriors(i,j).num_events,0:8,'Normalization','countdensity')
            axis off
% 
%             ylim([0 40]*10)
%             xlim([0 8])
% %             
%             figure(9993)
%             subplot1((i-1)*size(spatial_posteriors,2) + j);
%             histogram(spatial_posteriors(i,j).tau1,[5:5:50],'Normalization','countdensity')
%             axis tight
%             axis off
%             ylim([0 10000])
% %             xlim([5 50])
%             
            figure(99944)
            subplot1((i-1)*size(spatial_posteriors,2) + j);
            histogram(spatial_posteriors(i,j).tau2,[30:5:600],'Normalization','countdensity')
            axis tight
            axis off
            ylim([0 40]*10*3)
            xlim([30 600])
            

            figure(9995)
            subplot1((i-1)*size(spatial_posteriors,2) + j);
            histogram(spatial_posteriors(i,j).times,(0:.001:.1)*20000,'Normalization','countdensity') %
            axis off
            ylim([0 25]*10*3)
            xlim([0 .1]*20000)
            
%         end
    end
end

figure(9991)
subplot1(1);
title('amps')
% 
% figure(9992)
% subplot1(1)
% title('num events')
% 
% figure(9993)
% subplot1(1)
% histogram(spatial_posteriors(i,j).tau1)
% title('tau on')
% 
figure(99944)
subplot1(1)
title('tau off')
% 
figure(9995)
subplot1(1)
title('times')

    
    
function spatial_posteriors = build_spatial_posteriors_bivariate(posteriors_grid, neighborhood_size)

spatial_posteriors = struct();


figure(9991)
subplot1(size(posteriors_grid,1), size(posteriors_grid,2), 'Gap', [.005 .005], 'XTickL', 'Margin', 'YTickL', 'Margin');
% figure(9992)
% subplot1(size(posteriors_grid,1), size(posteriors_grid,2), 'Gap', [.005 .005], 'XTickL', 'Margin', 'YTickL', 'Margin');
% figure(9993)
% subplot1(size(posteriors_grid,1), size(posteriors_grid,2), 'Gap', [.005 .005], 'XTickL', 'Margin', 'YTickL', 'Margin');
% figure(99944)
% subplot1(size(posteriors_grid,1), size(posteriors_grid,2), 'Gap', [.005 .005], 'XTickL', 'Margin', 'YTickL', 'Margin');
% figure(9995)
% subplot1(size(posteriors_grid,1), size(posteriors_grid,2), 'Gap', [.005 .005], 'XTickL', 'Margin', 'YTickL', 'Margin');

for i = 1:size(posteriors_grid,1)
    for j = 1:size(posteriors_grid,2)
        
        if ~isempty(posteriors_grid{i,j})
            
            spatial_posteriors(i,j).amps = [];
            spatial_posteriors(i,j).num_events = [];
            spatial_posteriors(i,j).tau1 = [];
            spatial_posteriors(i,j).tau2 = [];
            spatial_posteriors(i,j).times = [];
            
            for ii = -neighborhood_size:neighborhood_size
                for jj = -neighborhood_size:neighborhood_size
                     if ~(i+ii < 1 || i+ii > size(posteriors_grid,1) || j+jj < 1 || j+jj > size(posteriors_grid,2))
                         
                         
                         spatial_posteriors(i,j).amps = [spatial_posteriors(i,j).amps [posteriors_grid{i+ii,j+jj}.amps]];
                         spatial_posteriors(i,j).num_events = [spatial_posteriors(i,j).num_events [posteriors_grid{i+ii,j+jj}.num_events]];
                         spatial_posteriors(i,j).tau1 = [spatial_posteriors(i,j).tau1 [posteriors_grid{i+ii,j+jj}.tau1]];
                         spatial_posteriors(i,j).tau2 = [spatial_posteriors(i,j).tau2 [posteriors_grid{i+ii,j+jj}.tau2]];
                         spatial_posteriors(i,j).times = [spatial_posteriors(i,j).times [posteriors_grid{i+ii,j+jj}.times]];
                         
                     end
                end
            end

            
            (i-1)*size(posteriors_grid,2) + j

                        
            figure(9991)
            subplot1((i-1)*size(posteriors_grid,2) + j);
            if ~isempty(spatial_posteriors(i,j).amps)
                hist3([spatial_posteriors(i,j).amps' spatial_posteriors(i,j).times'],'Edges',{[40:20:150],[0:.005:.1]*20000},'EdgeAlpha',0)
                colormap(hot)
                set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
            end
            view(3); 
            axis off
%             ylim([0 10000])
%             xlim([20 250])
            
%             figure(9992)
%             subplot1((i-1)*size(posteriors_grid,2) + j);
%             hist3([spatial_posteriors(i,j).num_events' spatial_posteriors(i,j).times'],'Edges',{0:8,[0:.005:.1]*20000})
%             axis off

%             ylim([0 10000])
%             xlim([0 7])
%             
%             figure(9993)
%             subplot1((i-1)*size(posteriors_grid,2) + j);
%             hist3([spatial_posteriors(i,j).tau1' spatial_posteriors(i,j).times'],'Edges',{[5:5:50],[0:.005:.1]*20000})
%             axis tight
%             axis off
% %             ylim([0 10000])
% %             xlim([5 50])
%             
%             figure(99944)
%             subplot1((i-1)*size(posteriors_grid,2) + j);
%             hist3([spatial_posteriors(i,j).tau2' spatial_posteriors(i,j).times'],'Edges',{[30:20:300],[0:.005:.1]*20000})
%             axis tight
%             axis off
%             ylim([0 10000])
%             xlim([30 200])
            

%             figure(9995)
%             subplot1((i-1)*size(posteriors_grid,2) + j);
%             histogram(spatial_posteriors(i,j).times,[0:.005:.1]*20000)
%             axis off
%             ylim([0 10000])
%             xlim([0 .1]*20000)
            
        end
    end
end

figure(9991)
subplot1(1);
title('amps')

% figure(9992)
% subplot1(1)
% title('num events')
% 
% figure(9993)
% subplot1(1)
% histogram(spatial_posteriors(i,j).tau1)
% title('tau on')
% 
% figure(99944)
% subplot1(1)
% title('tau off')

% figure(9995)
% subplot1(1)
% title('times')

    
    
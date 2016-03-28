function [feature_mat] = plot_event_features_grid(results_grid,feature_min)


feature_mat = [];

for i = 1:size(results_grid,1)
    for j = 1:size(results_grid,2)

            
            
            
            results = results_grid{i,j};
            
            for ii = 1:length(results)
                
                taus = zeros(0,2);
                good_events = [];
            
                for k = 1:length(results(ii).trials.tau{results(ii).map_ind})
                    if ~isempty(results(ii).trials.tau{results(ii).map_ind}{k})
                        if results(ii).trials.amp{results(ii).map_ind}(k) > 50

                            good_events = [good_events k];
                            taus = [taus; results(ii).trials.tau{results(ii).map_ind}{k}];

                        end
                    end

                end
                
                if ~isempty(good_events)
                    new_features = [results(ii).trials.amp{results(ii).map_ind}(good_events)' taus results(ii).trials.times{results(ii).map_ind}(good_events)'/20000 repmat([i j],length(good_events),1)];
                    feature_mat = [feature_mat; new_features];
                end
                
            
            end

    end
end

this_feature = feature_mat(:,4);
max_amp = max(this_feature) - feature_min;
max_amp = .03;
colors = [(this_feature - feature_min)/max_amp repmat([0 0],size(feature_mat,1),1)];

figure
scatter3(feature_mat(:,5),feature_mat(:,6),this_feature,[],colors,'filled')
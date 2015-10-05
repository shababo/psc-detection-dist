%% save parameter files for cluster

p_spike = [1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2];
a_min = [0 1e-2 1e-1 1e0 1e1 1e2 5e-2 5e-1 5e0 5e1 5e2];
num_sweeps = [1e1 1e2 1e3 1e4];

savefile_basename = 'data/cluster-param-files/params-pspike-%0.0e-amin-%0.0e-num_sweeps-%0.0e.mat';

for i = 1:length(p_spike)
    for j = 1:length(a_min)
        for k = 1:length(num_sweeps)
            params.cluster = 1;
            params.p_spike = p_spike(i);
            params.a_min = a_min(j);
            params.num_sweeps = num_sweeps(k);
            params = get_params(params);
            savefile = sprintf(savefile_basename,params.p_spike,params.a_min,params.num_sweeps);
            save(strrep(savefile,'+',''),'params')
            clear params
        end
    end
end
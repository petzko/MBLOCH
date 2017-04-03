clear

simfilenames = {'sim_scenarios/sec1.set','sim_scenarios/sec2.set','sim_scenarios/sec3.set',...
    'sim_scenarios/sec4.set','sim_scenarios/sec5.set','sim_scenarios/sec6.set'};
savefilenames = {'name1','name2','name3','name4','name5','name6'};
interpDataFile = 'fitted_data_OPTICA.mat';

Ni = length(simfilenames); 
workers = Ni;

spmd(workers)
    tidx = labindex;    
    acmodelocking(simfilenames{tidx},interpDataFile, savefilenames{tidx},false);
end

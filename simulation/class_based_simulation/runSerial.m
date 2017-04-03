clear

simfilename = {'sim_scenarios/sec1.set'};
savefilename = {'name1'};
interpDataFile = 'fitted_data_OPTICA.mat';

Ni = length(simfilename); 
acmodelocking(simfilename{1},interpDataFile, savefilename{1},true);


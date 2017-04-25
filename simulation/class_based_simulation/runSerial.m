close all;clear;

simfilename = {'sim_scenarios/sec4.set'};
savefilename = {'name4'};
interpDataFile = 'fitted_data_OPTICA.mat';

Ni = length(simfilename); 
acmodelocking(simfilename{1},interpDataFile,savefilename{1},true);


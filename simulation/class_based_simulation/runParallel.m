clear

simfilenames = {'sec1.set','sec2.set','sec3.set','sec4.set','sec5.set','sec6.set'};
savefilenames = {'name1','name2','name3','name4','name5','name6'};
interpDataFile = 'fitted_data_OPTICA.mat';

Ni = length(simfilenames); 
workers = Ni;

spmd(workers)
    tidx = labindex;    
    acmodelocking(simfilenames{tidx},interpDataFile, savefilenames{tidx});
end

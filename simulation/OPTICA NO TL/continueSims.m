clear; close all; clc;
in_folder = '.'; 
ratesfile = 'fitted_data.mat';

init = +1;

simfiles = {'OPTICAMAIN.sim'};
scenariofiles = {'12-1-verystrong.set','12-2-verystrong.set'};
executables = {@runsim_MB_TL_3lvl_FO};
workspacefiles  = {};

Ni = length(simfiles); 
Nj = length(scenariofiles); 
Nk =length(executables); 

workers = Ni*Nj*Nk;

spmd(workers)
    tidx = labindex;    
    i = floor((tidx-1)/((Nj)*(Nk)))+1; 
    j = floor(((tidx-1)-(i-1)*Nj*Nk)/Nk)+1;
    k = tidx-(j-1)*Nk-(i-1)*Nj*Nk;
    simfile = [in_folder '\' simfiles{i}];
    scenariofile = [in_folder '\' scenariofiles{j}];
    executable =  executables{k};
    savename = executable(scenariofile,simfile,'init',init,'workspace',workspacefiles{tidx});
end

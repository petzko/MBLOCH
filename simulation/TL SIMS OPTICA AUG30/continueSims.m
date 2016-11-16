clear; close all; clc;
in_folder = '.'; 
ratesfile = 'fitted_data.mat';

init = -1;

simfiles = {'sim_files/OPTICAMAIN.sim'};
scenariofiles = {'sim_files/12-4-0p1.set'};
executables = {@runsim_MB_TL_3lvl_FO};
workspacefiles  = {'CHCKPT_OPTICA(GHOSTCELL-newDRY-petzko)_(12-4-0p1)_N_TRANSMISSION_LINE_4000_FP'};

Ni = length(simfiles); 
Nj = length(scenariofiles); 
Nk =length(executables); 

workers = Ni*Nj*Nk;

% spmd(workers)
    tidx = labindex;    
    i = floor((tidx-1)/((Nj)*(Nk)))+1; 
    j = floor(((tidx-1)-(i-1)*Nj*Nk)/Nk)+1;
    k = tidx-(j-1)*Nk-(i-1)*Nj*Nk;
    simfile = [in_folder '\' simfiles{i}];
    scenariofile = [in_folder '\' scenariofiles{j}];
    executable =  executables{k};
    savename = executable(scenariofile,simfile,ratesfile,'init',init,'workspace',workspacefiles{tidx});
% end

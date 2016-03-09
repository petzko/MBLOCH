clear; close all; clc;
in_folder = '..\MB Inputs';  
init = -1;

simfiles = {'mlvl_11p0.sim'};
scenariofiles = {'standarddeph.set'};
executables = {@runsim_FP_RNFD_multilevelHTB_quadphase,@runsim_FP_RNFD_multilevelHTB_dispcompPHASE,@runsim_FP_RNFD_multilevelHTB};
workspacefiles  = {...
'qcl183s(MLVL,11kVpCm)_(standard-SHB)_WITH_QUADRATIC_PHASE_N_4000_FP_RT_1000.mat',...
'qcl183s(MLVL,11kVpCm)_(standard-SHB)_WITH_phase_compensation_after500_N_4000_FP_RT_1000.mat',...
'qcl183s(MLVL,11kVpCm)_(standard-SHB)_N_4000_FP_2000.mat'};

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

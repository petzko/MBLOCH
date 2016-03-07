clear; close all; clc;
in_folder = 'MB Inputs'; 

init = -1;

simfiles = {'mlvl_11p0.sim'};
scenariofiles = {'weakdeph.set','strongdeph.set','standarddeph.set'};
executables = {@runsim_FP_RNFD_multilevelHTB_quadphase,@runsim_FP_RNFD_multilevelHTB_dispcompperfectComp};
workspacefiles  = {...
'CHCKPT_qcl183s(MLVL,11kVpCm)_(weak-SHB)WITH-QUAD-PHASE_N_4000_FP.mat',...
'CHCKPT_qcl183s(MLVL,11kVpCm)_(weak-SHB)_WITH_PERFECT_PHASE_COMP_N_4000_FP.mat',...
'CHCKPT_qcl183s(MLVL,11kVpCm)_(strong-SHB)WITH-QUAD-PHASE_N_4000_FP.mat',...
'CHCKPT_qcl183s(MLVL,11kVpCm)_(strong-SHB)_WITH_PERFECT_PHASE_COMP_N_4000_FP.mat',...
'CHCKPT_qcl183s(MLVL,11kVpCm)_(standard-SHB)WITH-QUAD-PHASE_N_4000_FP.mat',...
'CHCKPT_qcl183s(MLVL,11kVpCm)_(standard-SHB)_WITH_PERFECT_PHASE_COMP_N_4000_FP.mat'};

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
    executable = executables{k};
    savename = executable(scenariofile,simfile,'init',init,'workspace',workspacefiles{tidx});
end

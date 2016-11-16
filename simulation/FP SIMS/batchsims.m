clear; close all; clc;
in_folder = '.';

init = +1;

simfiles = {'FL183S.sim'};
% scenariofiles = {'weakdeph.set','standarddeph.set','strongdeph.set'};
% executables = {@runsim_FP_RNFD_multilevelHTB_quadphase,@runsim_FP_RNFD_multilevelHTB_dispcompperfectComp};
scenariofiles = {'FL183S.set'};
executables = {@runsim_FP_RNFD_multilevelHTB};

Ni = length(simfiles); 
Nj = length(scenariofiles); 
Nk = length(executables); 

workers = Ni*Nj*Nk;

% spmd(workers)
    tidx = labindex;    
    i = floor((tidx-1)/((Nj)*(Nk)))+1; 
    j = floor(((tidx-1)-(i-1)*Nj*Nk)/Nk)+1;
    k = tidx-(j-1)*Nk-(i-1)*Nj*Nk;
%     simfile = [in_folder '\' simfiles{i}];
%     scenariofile = [in_folder '\' scenariofiles{j}];
    executable = executables{k};
    savename = executable(scenariofiles{j},simfiles{i},'init',init);
% end

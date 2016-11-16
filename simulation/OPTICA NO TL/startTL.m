clear; close all; clc;
in_folder = '.';

init = +1;

simfiles = {'OPTICAMAIN.sim'};
scenariofiles = {'12-1-off.set'};

Ni = length(simfiles);
Nj = length(scenariofiles);
Nk = 1;

workers = Ni*Nj*Nk;

tidx = labindex;
i = floor((tidx-1)/((Nj)*(Nk)))+1;
j = floor(((tidx-1)-(i-1)*Nj*Nk)/Nk)+1;
k = tidx-(j-1)*Nk-(i-1)*Nj*Nk;
simfile = [in_folder '/' simfiles{i}];
scenariofile = [in_folder '/' scenariofiles{j}];
savename = runsim_MB_TL_3lvl_FO(scenariofile,simfile,'init',init);

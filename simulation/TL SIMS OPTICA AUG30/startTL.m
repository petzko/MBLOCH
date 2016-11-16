clear; close all; clc;
in_folder = '.'; 
ratesfile = 'fitted_data.mat';

init = +1;

simfiles = {'sim_files\OPTICAMAIN.sim'};
% scenariofiles = {'sim_files\12-1-0p04.set',...
%     'sim_files\12-3-verystrong.set','sim_files\12-1-off.set'};

scenariofiles = {'sim_files\12-4-0p1.set'};


Ni = length(simfiles); 
Nj = length(scenariofiles); 
Nk = 1; 

workers = Ni*Nj*Nk;
% spmd(workers)

for file_idx = 1:length(scenariofiles)
    tidx = labindex;    
    i = floor((tidx-1)/((Nj)*(Nk)))+1; 
    j = floor(((tidx-1)-(i-1)*Nj*Nk)/Nk)+1;
    k = tidx-(j-1)*Nk-(i-1)*Nj*Nk;
    simfile = [in_folder '/' simfiles{i}];
%     scenariofile = [in_folder '/' scenariofiles{j}];
    scenariofile = [in_folder '/' scenariofiles{file_idx}];
    savename = runsim_MB_TL_3lvl_FO(scenariofile,simfile,ratesfile,'init',init);
end
% end
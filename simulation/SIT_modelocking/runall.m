clear; clc; close all; 

simDataFiles = {'SITMODELOCKING1.sim',... 
                'SITMODELOCKING1p5.sim',...
                 'SITMODELOCKING2.sim','SITMODELOCKING3.sim',...
                 'SITMODELOCKING3p6.sim'}

simDataFiles = { 'SITMODELOCKING1.sim'}
             
for file_idx = 1:length(simDataFiles)
    simDataFile = simDataFiles{file_idx} 
    sit_modelocking;     
end
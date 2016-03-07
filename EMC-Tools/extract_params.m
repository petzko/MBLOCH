function [ E_v,W_v ] = extract_params( parent_folder,base_vals,names,levels,wfs)
% 
% Extracts the wavefunctions and the energies from monte carlo simulation results stored in 
% all folders within the given parent_folder. the base_vals is an
% array specifying the values of the fixed variable (e.g. voltage or temperature)
% that determines the simulation content within each subfolder. the names
% of the subfolders are specified in the cell array 'names'
% note that the length of base_vals has to be the same as the length of the
% names array

%%% for example if the base_vals is a grid specifying different voltages,
%%% say [5,5.5,6,6.5,7,7.5 ... ] then in that same order the names array
%%% should contain the strings for subfolders with MC simulation results
%%% for each corresponding bias. 

%%% The input "levels" array, is an integer array specifying which
%%% wavefunctions to extract the information about. It is implicitly
%%% assumed that the wavefunctions are ordered in increasing energy!
%%% the wfs array is used as input to the scan function and tells it which
%%% and how many wavefuncitons there are. for 5 wavefunctions, wfs = 1:5

%%% The method returns a matrix of size length(base_vals)xlength(levels) 
%%% with the energies, and a 3D array length(base_vals)x(length(levels)xlength(levels))
%%% with the corresponding interlevel scattering rates. 



assert(length(base_vals) == length(names),'Arrays dimension mismatch. Please check your input data and try agian. ');

% extract only the directories within the parent_folder

W_v = zeros(length(base_vals),length(levels),length(levels)); 
E_v = zeros(length(base_vals),length(levels)); 


for f = 1:length(base_vals)

    cd([parent_folder '\' names{f}]);
    curr_dir = pwd; 
    %get the energies.. and the rates
    Ens = load('E_BOUND'); 
    % sort them 
    [sortedEns, order]  = sort(Ens(1:end,2));
    sortedNrs = Ens(order,1);
    %sort the levels
    Energy_lvls = sortedNrs(levels);
    Energies = sortedEns(levels);
    
    E_v(f,:) = Energies'; % take the transponent! 
    for i=1:length(Energy_lvls)
        for j = 1:length(Energy_lvls)
            if j~=i
                global sf;
                scan(curr_dir,Energy_lvls(i), Energy_lvls(j), wfs,[],3);
                W_v(f,i,j) = sf;
            end
        end
    end
    cd([parent_folder]);
end






end


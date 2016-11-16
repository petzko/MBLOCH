classdef SECTION
    %SECTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % density matrix list of vectors -> it is initialized by dm_model
        rho_; 
        % density matrix model 
        dm_model_; 
        % Hamiltonian 
        H_; 
        % bias vector 
        v_; 
        % current vector
        i_; 
        W_; % scattering rates matrix 
        
    end
    
    methods
        function [] = set_scattering_matrix(obj,new_W)
            % this function updates the scattering rates matrix between the laser levels
        end
        
        % this function sets the density matrix model used for the current
        % section -> 
        function [] = set_dm_model(obj,dm_model)
            
        end
    end
    
end


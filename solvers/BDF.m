classdef BDF < handle
    % Implements a BACKWARD - DIFFERENTIATION - FORMULAs for ODE's
    % 
    
    properties
        coefs; %% coeficients
        coef_0; % coef in front of y^n;
        m; %% nr of steps... of the multistep BDF method 
        N; %% nr of points on the spatial grid..
        %data container storing the solution in the previous "m" steps;
        %data is a R^(Nxm) matrix ...
        data;
        % index that cyclicly iterates over the numbers 1,2 .... m to indicate which 
        % data storage location should be replaced after the next
        % computational step
        storage_idx = 1; 
        %counts how many iterations have been made so far... 
        iter_ctr = 1;
        
    end
    
    methods
 
        function obj = BDF(nr_steps,nr_pts,init_data)
            % do some error checking ... 
            assert(nr_steps >=1,'please set a positive nr of steps for the MULTIstep algorithm to run correctly.'); 
            assert(nr_pts >= 2,'please set a greater value for your grid size');
            [ n1,n2 ]  = size(init_data);
            string = ['multistep algorihtm initialization failure.' ... 
                'initial data size and specified grid size do not agree...'] ; 
            assert(n1 == nr_pts && n2 == 1, string);
           
            %if all data is correctly initalized, proceed with the
            %constructor... 
            obj.m = nr_steps;
            obj.N = nr_pts;
            obj.data = zeros(obj.N,obj.m);
            obj.data(:,1) = init_data; % this sets - up the inital data...
            obj.coefs = zeros(obj.m,1);
            
        end
        
        function res = make_step(obj,rhs,dt)  
            
            % Make a single propagation step of size "dt"  from tn = n*dt -> to tn+1 = (n+1)*dt, using the stored
            % previous data. input arguments are the solver object itself,
            % the right hand side of the equation evaluated at the current timestep (tn) and the timestep size dt
            
            [n1,n2] = size(rhs); 
            %check if vector dimensions are consice.. 
            err_msg = ['Cannot evolve equation ' ... 
            'Rhs vector dimension does not agree with solution vector dimension.' ];
            assert(n1 == obj.N,err_msg); 
            
           %%% this ensures that we get the right initial conditions    
            if(obj.iter_ctr < obj.m)            
            
               [tmp , obj.coef_0] = get_coeffs_multistep(obj.iter_ctr);
               obj.coefs = [ zeros(obj.m-obj.iter_ctr,1) ; tmp  ]; % padd the coefs with zeros... 
               
            elseif(obj.iter_ctr == obj.m)
               [obj.coefs, obj.coef_0] = get_coeffs_multistep(obj.iter_ctr);
            end
            
           % after we complete more than m iterations the program will not get into this branch anymore...  
            
            % calculate and return the result 
            res = dt*rhs / obj.coef_0 ;
            for k = 1:obj.m
                  %in the iter_ctr'th iteration... 
                  f_k_s = mod(obj.m - obj.iter_ctr + k-1,obj.m) + 1; % cycle from left to right... 
                  %to the data stored at position k! 
                  res = res - obj.coefs(f_k_s)*obj.data(:,k);
            end
            
            %store the currently computed result at place storage_idx... 
            obj.storage_idx = mod(obj.storage_idx,obj.m) + 1; %% storage_idx is a periodic index with period m
            obj.data(:,obj.storage_idx) = res;  % finally store the result... 

            % after the step is made, increment the internal iterations counter of the object...
            obj.iter_ctr = obj.iter_ctr + 1;
            
            
        end
        function solution = get_latest_solution(obj)
         
            solution = obj.data(:,obj.storage_idx);
           
           end
        
        function idx = getIdx(storage_ind, k,m,iter_ctr)
            
            if(iter_ctr <=m)
                idx = k; 
            else
                idx =  mod(storage_ind + k-1,m) + 1;
            end
                
                
        end
    end
    
end


function res = m_choose_i(m,i)
    
    res = 1; 
    for p = 1:i
       res = res * (m-p+1)/p ; 
    end

end

function [coefs,c_0] = get_coeffs_multistep(nr_steps)

            c_0 = 0; 
            coefs = zeros(nr_steps,1);
            for j = 1:nr_steps
                c_0 = c_0 + 1.0/j; 
            end
            
            %compute the rest of the coeffs... 
            pow = -1;
            m = nr_steps;
            for i = 1:m
                %store the coefs in reverse order..
                coefs(m-i+1) = 1/c_0 * pow*m_choose_i(m,i)/i;
                pow = pow*-1;
            end
   

end

classdef MS_O1 < handle
    % Implements an explicit lineear multistep algorihtm for ODE's
    
    
    properties
        
        % specify the coefficients for 1,2,3,4 and 5 step adams bashforth
        % method... so coefs(k,p) gives you the pth coefficient for a
        % k-step adams bashforth system...
        MAX_STEPS = 5
        %%% coeffs are stored in the following manner "zeros padding <---> oldest <---> newest "
        coefs; % = [ 0 0 0 0 1; 0 0 0 -1/2 3/2; 0  0 5/12 -4/3 23/12; 0 -3/8 37/24 -59/24 55/24 ; 251/720 -637/360 109/30 -1387/360 1901/720];
        coef_0 = 1; % coef in front of y^n;
        m; %% nr of steps...
        N; %% nr of points on the spatial grid..
        %data container storing the rhs of the previous "m" steps;
        %data is a R^(Nxm) matrix ...
        data; %
        % a vector of size N which stores the solution at the previous time
        % step...
        prev_solution;
        %counts how many iterations have been made so far...
        iter_ctr = 1;
        storage_idx; 
        retrieval_order; 
        
    end
    methods
        function obj = MS_O1(nr_steps,nr_pts,init_rhs_data,init_solution)
            
            % do some error checking ...
            assert(nr_steps >=1,'please set a positive nr of steps for the MULTIstep algorithm to run correctly.');
            assert(nr_steps <= 5,['unfortunately current implementation does not support more then a 5 step adams bashforth algorithm!'...
                ' Please specify a new nr of steps parameter and try again,']);
            
            [rowdim,coldim] = size(init_rhs_data);
            % do some assertions
            assert((rowdim >=0 && rowdim <=(nr_steps-1)), 'Incorrectly specified initial data! Please try again');
            
            %if all data is correctly initalized, proceed with the
            %constructor...
            
            
            
            obj.m = nr_steps;
            obj.N = nr_pts;
            %
            obj.data = zeros(obj.m,obj.N);
            for k = 1:rowdim
                obj.data(k+1,:) = init_rhs_data(rowdim - k+1,:);
            end
            
            % solution at step n2
            obj.coefs = get_coeffs(obj.m);
            obj.prev_solution = init_solution; % this sets - up the inital data...
            obj.iter_ctr = rowdim+1; 
            
        end
        
        function res = make_step(obj,prev_rhs,dt)
         
           % Make a single propagation step of size "dt"  from tn = n*dt -> to tn+1 = (n+1)*dt, using the stored
            % previous data. input arguments are the solver object itself,
            % the right hand side of the equation evaluated at the current timestep (tn) and the timestep size dt
            
            
            % replace the oldest rhs with the newest rhs
            obj.data(1,:) = prev_rhs;  % finally store the result...
            
            % if we are still filling up the obj.data with initial data...
            if(obj.iter_ctr <= obj.m)
                tmp =  get_coeffs(obj.iter_ctr);
                obj.coefs = [tmp ;  zeros(obj.m-obj.iter_ctr,1) ]; % padd the coefs with zeros...
            end
            
            % now compute the solution
            res =  obj.prev_solution;
            for k = 1:obj.m
                res = res + dt*obj.coefs(k)*obj.data(k,:);
            end
            
            % after the step is made, increment the internal iterations counter of the object...
            obj.iter_ctr = obj.iter_ctr + 1;
            obj.prev_solution = res; % save solution y_n+1 as the previous solution

            %circularly shift the matrix so that the oldest rhs comes on
            %first spot.
            obj.data = circshift(obj.data,1);
        end
        
        function solution = get_latest_solution(obj)
            solution  = obj.prev_solution;
        end
    end
    
end

function coefs =  get_coeffs(step)
if(step == 1)
    coefs =  1;
elseif(step == 2)
    coefs = [3/2; -1/2 ];
elseif(step == 3)
    coefs =  [23/12; -4/3 ; 5/12 ];
elseif(step == 4)
    coefs = [55/24;-59/24;37/24; -3/8 ];
elseif(step == 5)
    coefs = [ 1901/720; -1387/360; 109/30; -637/360; 251/720 ];
end


end

classdef RNFDSolver < handle
    % Finite Difference wave eqn solver based on the paper of Risken &
    % Nummedal from 1968 th or sth... 
    
    properties        
        PROPERTIES; % a small GPU ARRAY containing Nr pts, Direction, Grid_SPACING and Velocity properties
        SOLUTION; % this is the solution vector 
    end
    
    methods
        function obj = RNFDSolver(N,dx,direction,velocity, U_0)
       
            assert(N >= 2, ['You have specified too few grid points. Please ' ...
            'check the input aparameters and try again.']); 
            assert(direction ~= 0, ['Please specify valid wave propagation direction parameter. > 0 '...
                'for forward and < 0 for backward propagating wave.']);
            
            [n1,n2] = size(U_0); 
            n1 = max(n1,n2);
            assert(N == n1, 'Initial data size and specified grid size do not agree. please check your input params and try again..');
            obj.PROPERTIES = zeros(4,1);
            obj.PROPERTIES(1) = N; 
            obj.PROPERTIES(2) = dx; 
            obj.PROPERTIES(3) = abs(velocity);
            obj.PROPERTIES(4) = direction / abs(direction);
            obj.SOLUTION = U_0;
            % send the props and the solution arrays to reside on the GPU
            obj.PROPERTIES = gpuArray(obj.PROPERTIES);
            obj.SOLUTION = gpuArray(U_0);
        end
        function solution = get_latest_solution(obj)
            solution = obj.SOLUTION;
        end

        
        function res = make_step(obj,F,F_t,k,dt)
         
                %%% The solver assumes an equation of the form 
                %%% D_t(E) = -/+ D_x(E) + f(x,t) + k E 
                %%% 
                
                
                E = obj.SOLUTION; % alias for the solution vector
                N = obj.PROPERTIES(1); % alias for the grid size
                dx = obj.PROPERTIES(2);
                c = obj.PROPERTIES(3);
                assert(c*dt == dx, 'RN FDsolver make step failed. Magic time step not selected.'); 
               
                dt_2 = (dt^2)/2;
                E_x = (E(2:N) - E(1:N-1))/dx;
                F_x = (F(2:N) - F(1:N-1))/dx;
                
                if(obj.PROPERTIES(4) > 0)
                     obj.SOLUTION(2:N) = E(1:N-1) +dt_2*( F_t(2:N) -c*F_x - 2*k*c*E_x + k*F(2:N)+k^2*E(2:N)) +dt*(F(2:N)+k*E(2:N)); 
                else    
                     obj.SOLUTION(1:N-1) = E(2:N) +dt_2*( F_t(1:N-1) + c*F_x + 2*k*c*E_x + k*F(1:N-1)+k^2*E(1:N-1)) +dt*(F(1:N-1)+k*E(1:N-1)); 
                end
                
                res = obj.SOLUTION;
                
        end
        
        function sol = set_bdry(obj, Lvalue,Rvalue)
             
            if(~strcmp('no',Lvalue))
               obj.SOLUTION(1) = Lvalue;
            end
            
            if(~strcmp('no',Rvalue))
                obj.SOLUTION(end)= Rvalue;
            end
            
            sol = obj.SOLUTION; % return the updated solution... 
        end
            
            
    end
    
end


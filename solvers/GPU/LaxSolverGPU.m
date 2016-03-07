classdef LaxSolverGPU < handle
    %SPLITSTEPSOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        PROPERTIES; % a small GPU ARRAY containing Nr pts, Direction, Grid_SPACING and Velocity properties
        SOLUTION; % this is the solution vector 
       
    end
    
    methods
        function obj = LaxSolverGPU(N,dx,direction,velocity, U_0)
       
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
            obj.PROPERTIES(3) = velocity;
            obj.PROPERTIES(4) = direction / abs(direction);
            
            % upload the props and the initial data on GPU memory
            obj.PROPERTIES = gpuArray(obj.PROPERTIES);
            obj.SOLUTION = gpuArray(U_0);
        end
        
        function solution = get_latest_solution(obj)
            solution = obj.SOLUTION;
        end
        
        % for one lax wendroff step one needs to have the rhs vector
        % at the current step as well as the time and space derivatives of
        % the rhs vector... NOTE THAT THIS SOLVER DOES NOT IMPLEMENT ANY
        % BOUNDARY CONDITIONS!!!!!!!!! 
        function res = make_step(obj,f,f_t,k,dt)
           
            % THE SOLVER is adapted to solving an equation in the form: 
            %   D_t E(x,t) = -/+ c D_x E(x,t) + f(x,t) + k E
            if(obj.PROPERTIES(4) > 0 ) 
                U= obj.SOLUTION; % alias for the solution vector
                N = obj.PROPERTIES(1); % alias for the grid size
                c = obj.PROPERTIES(3); % alias for the phase velocity
                dx = obj.PROPERTIES(2); % alias for the grid size
                f_x = (f(3:N) - f(1:N-2))/(2*dx);
                Ux = (U(3:N) - U(1:N-2))/(2*dx); % gpu array already...
                Uxx = (U(3:N) -2* U(2:N-1) + U(1:N-2))/(dx*dx);
                obj.SOLUTION(2:N-1) = U(2:N-1) + dt*(-c*Ux+f(2:N-1)+k*U(2:N-1))+...
                    (dt*dt/2)*((c^2)*Uxx+ (f_t(2:N-1)-c*f_x) +k*(-2*c*Ux +f(2:N-1)+k*U(2:N-1)));
            else
                V = obj.SOLUTION; %alias for the old solution !!!
                N = obj.PROPERTIES(1); % alias for the grid size
                c = obj.PROPERTIES(3); % alias for the phase velocity
                dx = obj.PROPERTIES(2); % alias for the grid size
                f_x = (f(3:N) - f(1:N-2))/(2*dx);
                Vx = (V(3:N) - V(1:N-2))/(2*dx);
                Vxx = (V(3:N) -2* V(2:N-1) + V(1:N-2))/(dx*dx);
                obj.SOLUTION(2:N-1,1) = V(2:N-1,1) + dt*(c*Vx+f(2:N-1)+k*V(2:N-1)) + ...
                (dt*dt/2)*((c^2)*Vxx + (f_t(2:N-1)+c*f_x) +k*(2*c*Vx  + f(2:N-1) + k*V(2:N-1)));
            end
            
             res = obj.SOLUTION;
        
        end
        function sol = set_bdry(obj, Lvalue , Rvalue )

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


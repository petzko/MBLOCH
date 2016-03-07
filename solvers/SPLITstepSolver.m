classdef SplitstepSolver < handle
    %SPLITSTEPSOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        GRID_SPACING; % grid spacing size 
        VELOCITY; % phase velocity of the wave under investigation 
        SOLUTION; % this is the solution vector 
        NR_PTS; % specifies the nr of points on the grid...
        LOSSES; % specify the loss coefficient
        DIR; % +1 if wave is in positive direction and -1 if wave is in negative direction ... 
        PHASE_SHIFT;
        SCHEME;
    end
    
    methods
        function obj = SplitstepSolver(Nr_pts,dx,direction,loss,velocity, U_0)
       
            assert(Nr_pts >= 2, ['You have specified too few grid points. Please ' ...
            'check the input aparameters and try again.']); 
            assert(direction ~= 0, ['Please specify valid wave propagation direction parameter. > 0 '...
                'for forward and < 0 for backward propagating wave.']);
            
            [n1,n2] = size(U_0); 
            n1 = max(n1,n2);
            assert(Nr_pts == n1, 'Initial data size and specified grid size do not agree. please check your input params and try again..');
            obj.NR_PTS = Nr_pts; 
            obj.GRID_SPACING = dx; 
            obj.SOLUTION = zeros(obj.NR_PTS,1);
            obj.SOLUTION = U_0;
            obj.DIR = direction / abs(direction); 
            obj.LOSSES = loss; % setup teh loss coefficient...            
            obj.VELOCITY = velocity; 
   
         
          
        end
        function solution = get_latest_solution(obj)
            solution = obj.SOLUTION;
        end
        % for one lax wendroff step one needs to have the rhs vector
        % at the current step as well as the time and space derivatives of
        % the rhs vector... NOTE THAT THIS SOLVER DOES NOT IMPLEMENT ANY
        % BOUNDARY CONDITIONS!!!!!!!!! 
        function res = make_step(obj,f,dt)
                
            %obtain aliases for linear losses, phase velocity and previous
            %solution
            
                l = obj.LOSSES;
                c = obj.VELOCITY;
                U = obj.SOLUTION;
                direction = obj.DIR;
                
                alpha = -direction*c*dt/obj.GRID_SPACING;
                phi = 2*pi*alpha/obj.NR_PTS;
                k_v = linspace(0,obj.NR_PTS-1,obj.NR_PTS)';
                obj.PHASE_SHIFT = exp(1i*phi*k_v);
                
            % do a symmetric strang splitting ...     
                U = ifft(fft(U).*obj.PHASE_SHIFT); % this is the analytical solution... 
%                 U = U + 0.5*dt*c*f -l*c*dt*U;
%                 U = ifft(fft(U).*obj.PHASE_SHIFT);
                obj.SOLUTION = U;
                res = obj.SOLUTION;
        
        end
        
        
        function sol = set_bdry(obj, LValue , RValue )
                obj.SOLUTION(1) = LValue;
                obj.SOLUTION(end) = RValue;
                sol = obj.SOLUTION; % return the updated solution... 
        end
     
            
            
    end
    
end


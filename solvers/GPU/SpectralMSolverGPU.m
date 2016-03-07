classdef SpectralMSolverGPU < handle
    %SPECTRALMSOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        PROPERTIES;
        N; % specifies the nr of points on the grid...
        DIR; % direction
        VELOCITY; % phase velocity of the wave under investigation 
        
        SOLUTION; % this is the solution vector 
        PROPAGATOR;
        M; %differentiation Mtx
    end
    
    methods
        function obj = SpectralMSolverGPU(N,direction,a,b,nr_steps,velocity,U_0)
              assert(N >= 2, ['You have specified too few grid points. Please ' ...
            'check the input aparameters and try again.']); 
            assert(direction ~= 0, ['Please specify valid wave propagation direction parameter. > 0 '...
                'for forward and < 0 for backward propagating wave.']);

            [n1,n2] = size(U_0); 
            assert(N == n1 && 1 == n2, ['Initial data dimensions incorrect. Please setup your solution vectors as column vectors, i.e.' ...
             'in Nx1 dimensional form. '] );
            assert(b > a , 'Spectral Solver initialization FAILED. Erroneous domain boundaries [a,b]!');
         
            obj.PROPERTIES = zeros(4,1);
            obj.PROPERTIES(1) = N; 
            obj.PROPERTIES(2) = nr_steps;
            obj.PROPERTIES(2) = velocity;
            obj.PROPERTIES(3) = direction / abs(direction);

            % upload the props and the initial data on GPU memory
            obj.PROPERTIES = gpuArray(obj.PROPERTIES);
            obj.SOLUTION = gpuArray(U_0);
            obj.PROPAGATOR = MSqueueGPU(nr_steps,N,[],U_0);

            n = N-1;
            obj.M = 2/(b - a)*obj.PROPERTIES(3)*obj.PROPERTIES(2)*obj.cheb(n);
            obj.M = gpuArray(obj.M);
                     
        end
        
          function res = make_step(obj,P,dt)
                P = gpuArray(P);
                rhs = obj.M*obj.SOLUTION + P;
                obj.SOLUTION = obj.PROPAGATOR.make_step(rhs,dt);
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
       
       function [D,x] = cheb(obj,N)
           
          if N==0, D=0; x=1; return, end
          x = cos(pi*(0:N)/N)'; 
          c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
          X = repmat(x,1,N+1);
          dX = X-X';                  
          D  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
          D  = D - diag(sum(D'));                 % diagonal entries
          
        end
    end
    
end


classdef HighOrderSolverImplicit < handle
    %SPLITSTEPSOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        DX; % grid spacing size 
        DT; % timestep size...
        VELOCITY; % phase velocity of the wave under investigation 
        SOLUTION; % this is the solution vector 
        PROPAGATOR
        NR_PTS; % specifies the nr of points on the grid...
        DIR; % +1 if wave is in positive direction and -1 if wave is in negative direction ... 
        C; % multiplication constant!
        M;
    end
    
    methods
      function obj = HighOrderSolverImplicit(N,dx,dt,direction,velocity, U_0)
       
            assert(N >= 2, ['You have specified too few grid points. Please ' ...
            'check the input aparameters and try again.']); 
            assert(direction ~= 0, ['Please specify valid wave propagation direction parameter. > 0 '...
                'for forward and < 0 for backward propagating wave.']);

            [n1,n2] = size(U_0); 
            assert(N == n1 && 1 == n2, ['Initial data dimensions incorrect. Please setup your solution vectors as column vectors, i.e.' ...
             'in Nx1 dimensional form. '] );
         
            obj.NR_PTS = N; 
            obj.DX = dx; 
            obj.DT = dt; 
            obj.SOLUTION = U_0;
            obj.DIR = direction / abs(direction); 
            obj.VELOCITY = velocity; 
            obj.C = -obj.DIR*obj.VELOCITY/(obj.DX);
%             n = N-1;
            n = N;
            obj.M = zeros(n);
            nr_steps =3;
            a = obj.C;
            M1 = zeros(n); M2 = zeros(n);
            
            M1(1,1) = 6; M1(1,2) = 18;
            for j = 2:n-1
                M1(j,j-1) = 1;
                M1(j,j) = 4;
                M1(j,j+1) = 1;
            end
            M1(n,n-1) = 18; M1(n,n)  = 6;

            M2(1,1) = -17; M2(1,2) = 9; M2(1,3) = 9; M2(1,4) = -1;
            for j = 2:n-1
                M2(j,j-1) = -3; M2(j,j) = 0; M2(j,j+1) = 3;
            end
            M2(n,n-3) = 1; M2(n,n-2) = -9; M2(n,n-1) = -9; M2(n,n) =17;
            obj.M = a*inv(M1)*M2;
            obj.PROPAGATOR = MSqueue(nr_steps,n,[],U_0);
                
    
      end
        
        function solution = get_latest_solution(obj)
            solution = obj.SOLUTION;
        end

        
        function res = make_step(obj,P,dt)
               
                V = obj.SOLUTION; % alias for the solution vector
                              
%                 if (obj.DIR > 0)
%                     obj.SOLUTION(1) = 0;
%                     rhs = obj.M*V(2:end,1) +obj.B*V(1,1) + P(2:end,1);
%                     obj.SOLUTION(2:end) = obj.PROPAGATOR.make_step(rhs,dt);
%                 else
%                     
%                     obj.SOLUTION(end) = 0;
%                     rhs = obj.M*V(1:end-1,1) +obj.B*V(end,1) + P(1:end-1,1);
%                     obj.SOLUTION(1:end-1,1) = obj.PROPAGATOR.make_step(rhs,dt);
%                     
%                 end

                   
                rhs = obj.M*V + P;
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
    end
    
end


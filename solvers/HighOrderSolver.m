classdef HighOrderSolver < handle
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
        M1;
        M2;
        B1;
        B;
    end
    
    methods
      function obj = HighOrderSolver(N,dx,dt,direction,velocity, U_0)
       
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
            obj.C = -obj.DIR*obj.VELOCITY/(12*obj.DX);
%             n = N-1;
            n = N;
            obj.M = zeros(n);
            nr_steps =3;
            a = obj.C;
            
            obj.M(1,1) = -25*a; obj.M(1,2) = 48*a ; obj.M(1,3) = -36*a; obj.M(1,4) = 16*a; obj.M(1,5) = -3*a;  
            obj.M(2,1) = -3*a ; obj.M(2,2) =  -10*a ; obj.M(2,3) =  18*a; obj.M(2,4) = -6*a;    obj.M(2,5) =1*a;

            for j = 3: n-2
                obj.M(j,j-2) = 1*a; 
                obj.M(j,j-1) = -8*a;
                obj.M(j,j) = 0*a;
                obj.M(j,j+1) = 8*a;
                obj.M(j,j+2) = -1*a;
            end

            obj.M(n-1,n-4) = -1*a; obj.M(n-1,n-3) = +6*a ; obj.M(n-1,n-2) =  -18*a ; obj.M(n-1,n-1) =  10*a; obj.M(n-1,n) = 3*a;
            obj.M(n,n-4) = 3*a ; obj.M(n,n-3) = -16*a ; obj.M(n,n-2) =  36*a ; obj.M(n,n-1) = -48*a; obj.M(n,n) = 25*a;
                
            obj.PROPAGATOR = MSqueue(nr_steps,n,[],U_0);
                
    
      end
        
        function solution = get_latest_solution(obj)
            solution = obj.SOLUTION;
        end

        
        function res = make_step(obj,P,dt)
                V = obj.SOLUTION; % alias for the solution vector
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


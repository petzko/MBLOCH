
classdef SpectralMSolver < handle
    %SPECTRALMSOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        DIR; % direction
        a,b; %boundaries of the domain...
        VELOCITY; % phase velocity of the wave under investigation
        SOLUTION; % this is the solution vector
        PROPAGATOR;
        N; % specifies the nr of points on the grid...
        M; %differentiation Mtx
        F;
        ode_options;
    end
    
    methods
        function obj = SpectralMSolver(N,direction,a,b,nr_steps,velocity,init_rhs_data,U_0)
            
            assert(N >= 2, ['You have specified too few grid points. Please ' ...
                'check the input aparameters and try again.']);
            assert(direction ~= 0, ['Please specify valid wave propagation direction parameter. > 0 '...
                'for forward and < 0 for backward propagating wave.']);
            
            [n1,n2] = size(U_0);
            assert(N == n1 && 1 == n2, ['Initial data dimensions incorrect. Please setup your solution vectors as column vectors, i.e.' ...
                'in Nx1 dimensional form. '] );
            
            obj.N = N;
            obj.SOLUTION = U_0;
            obj.DIR = direction / abs(direction);
            obj.VELOCITY = velocity;
            obj.a = a;
            obj.b = b;
            n = N-1;
            obj.M = -obj.DIR*obj.VELOCITY*2/(obj.b - obj.a)*SpectralMSolver.cheb(n);
            obj.PROPAGATOR = MS(nr_steps,N,init_rhs_data,U_0);
            
        end
        
        function res = make_step(obj,P,dt)
            
            
                rhs = obj.M*obj.SOLUTION + P;
                obj.SOLUTION = obj.PROPAGATOR.make_step(rhs,dt);
                res = obj.SOLUTION;
            
%             obj.F = P;
%             [TIMES, res] = ode45(@obj.diffeqn,[0 ,dt],obj.SOLUTION);
%             res = res(end,:)';
%             obj.SOLUTION = res;
                        
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
    methods (Static)
        function [D,x] = cheb(N)
            
            if N==0, D=0; x=1; return, end
            x = -cos(pi*(0:N)/N)';
            c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
            X = repmat(x,1,N+1);
            dX = X-X';
            D  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
            D  = D - diag(sum(D'));                 % diagonal entries
            
        end
        
        function x_prime = diffeqn(obj,t,x)
            x_prime = obj.M*x + obj.F;
        end
    end
    
end




classdef SpectralMSolverGPUKernel < handle
    %SPECTRALMSOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        MAX_STEPS = 5;
        MAX_THREADS_BLOCK = 1024;
    end
    
    properties
        
        N; % specifies the nr of points on the grid...
        m; % nr of steps for the multistep method..
        
        KERNEL; % class wrapper of the CUDA kernel that will perform the computation
        SOLUTION; % this is the solution vector
        D_N; %differentiation Mtx
        
        % specify the coefficients for 1,2,3,4 or 5 step adams bashforth
        
        %%% coeffs are stored in the following manner "zeros padding <---> oldest <---> newest "
        coefs; % = [ 0 0 0 0 1; 0 0 0 -1/2 3/2; 0  0 5/12 -4/3 23/12; 0 -3/8 37/24 -59/24 55/24 ; 251/720 -637/360 109/30 -1387/360 1901/720];
        
        
        %data container storing the rhs of the previous "m" steps;
        %F is a R^(Nxm) matrix ...
        F; %
        % a vector of size N which stores the solution at the previous time
        %counts how many iterations have been made so far...
        
        iter_ctr = 1;        
        storage_idx;
        order_idx;
        
        
    end
    
    methods
        function obj = SpectralMSolverGPUKernel(N,direction,a,b,nr_steps,velocity,U_0)
            
            
            % do some error checking ...
            assert(nr_steps >=1,'please set a positive nr of steps for the MULTIstep algorithm to run correctly.');
            assert(nr_steps <= 5,['unfortunately current implementation does not support more then a 5 step adams bashforth algorithm!'...
                ' Please specify a new nr of steps parameter and try again,']);
            
            assert(N >= 2, ['You have specified too few grid points. Please ' ...
                'check the input aparameters and try again.']);
            assert(direction ~= 0, ['Please specify valid wave propagation direction parameter. > 0 '...
                'for forward and < 0 for backward propagating wave.']);
            
            [n1,n2] = size(U_0);
            assert(N == n1 && 1 == n2, ['Initial data dimensions incorrect. Please setup your solution vectors as column vectors, i.e.' ...
                'in Nx1 dimensional form. '] );
            assert(b > a , 'Spectral Solver initialization FAILED. Erroneous domain boundaries [a,b]!');
            
            obj.N = N;
            obj.m = nr_steps;
            % upload the props and the initial data on GPU memory
            obj.SOLUTION = gpuArray(complex(U_0));
            n = N-1;
            
            obj.D_N = 2/(b - a)*direction/abs(direction)*velocity*obj.cheb(n);
            obj.D_N = gpuArray(reshape(obj.D_N,N*N,1));
            
            obj.KERNEL = parallel.gpu.CUDAKernel('SM_make_step.ptx','SM_make_step.cu');
            blockDim = [1024 1 1 ];
            gridDim = [ceil(N/blockDim(1)) 1  1];
            
            obj.KERNEL.GridSize = gridDim;
            obj.KERNEL.ThreadBlockSize = blockDim;            
            
            % propagator related stuff
            obj.F = complex(zeros(obj.N*obj.m,1));
            
            % solution at step n2
            obj.coefs = get_coeffs(obj.m);
            obj.iter_ctr = 1;
            tmp = 0:(obj.m-1);
            obj.order_idx =  gpuArray(int32(tmp));
            obj.storage_idx = int32(tmp);
            
        end
        
        function res = make_step(obj,P,dt)
      
            % if we are still filling up the obj.data with initial data...
            if(obj.iter_ctr <= obj.m)
                tmp =  get_coeffs(obj.iter_ctr);
                obj.coefs = gpuArray([tmp ;  zeros(obj.m-obj.iter_ctr,1) ]); % padd the coefs with zeros...
            end
            
            store_at = obj.storage_idx(1);        
            
%             __global__ void make_step(const double* D_N, const double2 * V, const double2* P,
% 							const double* c,
% 							const double dt, const int N, const int m, const int store_idx, 
% 							const int * order, 
% 							double2* F , double2* V_new)
            
            [obj.F,obj.SOLUTION] = feval(obj.KERNEL,obj.D_N,obj.SOLUTION,P,obj.coefs,...
                                        dt,obj.N,obj.m,store_at, ...
                                        obj.order_idx,obj.F,obj.SOLUTION); 
            
            res = obj.SOLUTION;
            
            % after the step is made, increment the internal iterations counter of the object...
            obj.iter_ctr = obj.iter_ctr + 1;
            obj.storage_idx = mod(obj.storage_idx -1,obj.m);
            obj.order_idx = mod(obj.order_idx+1,obj.m);
            
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


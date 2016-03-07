classdef MSqueueGPU < handle
    % Implements an explicit lineear multistep algorihtm for ODE's
    
    
    properties
        
        % specify the coefficients for 1,2,3,4 and 5 step adams bashforth
        % method... so coefs(k,p) gives you the pth coefficient for a
        % k-step adams bashforth system...
        MAX_STEPS = 5
        %%% coeffs are stored in the following manner "zeros padding <--->
        %%% oldest <---> newest " store the coefficients on GPU memory
        coefs; % = [ 0 0 0 0 1; 0 0 0 -1/2 3/2; 0  0 5/12 -4/3 23/12; 0 -3/8 37/24 -59/24 55/24 ; 251/720 -637/360 109/30 -1387/360 1901/720];
        coef_0 = 1; % coef in front of y^n;
        m; %% nr of steps...
        N; %% nr of points on the spatial grid..
        %data container storing the rhs of the previous "m" steps;
        dataQueue; % all data stored resides in GPU MEMORY
        
        % a vector of size N which stores the solution at the previous time
        % step...
        prev_solution; % on GPU memory
        %counts how many iterations have been made so far...
        iter_ctr;  % a counter specifying at which iteration/timestep is the solver so far!!!
        % if iter_ctr = n, this means that the NEXT call of make_step will be computing the solution at
        % the n-th timestep , i.e. tn = n*dt;
        
    end
    methods
        
        %  nr_steps:  the number of steps the algorithm will have ! Must be a
        %  number between 1 and 5
        %  nr_pts:  the number of data points in the solution vector
        %  initi_rhs_data: a matrix of dimensions nr_pts * n2 where n2 is the
        %       number of previous rhs evaluations provided!!! n2 has to be in the
        %       interval [1;5] 
        %  init_solution: a single data vector containing the solution of
        %  the latest!!!! time-step !!!
        % 
        %       !!! WARNING , for a obj.m steps adamsBashforth we
        %       need to have maximum of (obj.m-1) initial datas!!!
        %       for example, consider the 3 step adams-bashforth on the ODE
        %       dy/dt = f(t,y) 
        %       the initial rhs data needed is [f(0,y(0)) , f(dt,y(dt))] !!!
        %       then the initial solution needed in order to compute
        %       y(3*dt) is y(2*dt)!!!
        %       according to the formula, then the solution at time 3*dt
        %       will be given by
        %           y(3*dt) = y(2*dt) + dt*(f(2*dt,y(2*dt))*c1 + f(dt,y(dt))*c2+ f(0,y(0))*c3)
        %       , where f(2*dt,y(2*dt)) is an input argument to
        %       obj.make_step()
        %
        
        
        
        
        
        %       !!!ANOTHER WARNING - if Multiple initial conditions are supplied ny init_data, they have to be
        %       ordered from oldest to newest, i.e. if init_data(:,i),init_data(:,j) is the
        %       initial condition at time t_i and t_j respectively, and if i<j
        %       then  t_i < t_j !!!
        
        function obj = MSqueueGPU(nr_steps,nr_pts,init_rhs_data,init_solution)
            
            % do some error checking ...
            assert(nr_steps >=1,'please set a positive nr of steps for the MULTIstep algorithm to run correctly.');
            assert(nr_steps <= 5,['unfortunately current implementation does not support more then a 5 step adams bashforth algorithm!'...
                ' Please specify a new nr of steps parameter and try again,']);
        
            [ n1,n2 ]  = size(init_rhs_data);
            % do some assertions
            assert((n2 >=0 && n2 <=(nr_steps-1)), 'Incorrectly specified initial data! Please try again');
            
            %if all data is correctly initalized, proceed with the
            %constructor...
            
            obj.m = nr_steps;
            obj.N = nr_pts;
            obj.dataQueue = QueueGPU();
            
            % enqueue the oldest data points as zeros
            for k = 1:obj.m-n2
                obj.dataQueue.enqueue(zeros(nr_pts,1));
            end
            
            %enqueue the initial data to the begining of the queue...
            for k = 1:n2 
                obj.dataQueue.enqueue(init_rhs_data(:,k));
            end
            
            % solution at step n2
            
            obj.coefs = gpuArray(get_coeffs(obj.m));
            obj.prev_solution = gpuArray(init_solution); % this sets - up the inital data...
            obj.iter_ctr = n2+1; 
        end
        
        function res = make_step(obj,prev_rhs,dt)
     
            % Make a single propagation step of size "dt"  from tn = n*dt -> to tn+1 = (n+1)*dt, using the stored
            % previous data. input arguments are the solver object itself,
            % the right hand side of the equation evaluated at the current timestep (tn) and the timestep size dt
            
            %add the new rhs vector and remove the oldest one!!!
            obj.dataQueue.enqueue(prev_rhs);     
            obj.dataQueue.dequeue();
     
            % now compute the solution
            if (obj.iter_ctr <= obj.m)
                tmp = get_coeffs(obj.iter_ctr);
                obj.coefs = gpuArray([tmp;zeros(obj.m - obj.iter_ctr,1)]);
            end
                
            res =  obj.prev_solution;
            node = obj.dataQueue.getBackNode();
            k = 1;
            while (node ~= -1)
                data_k = node.getData();
                coef_k = obj.coefs(k);
                res = res + dt*coef_k*data_k;                
                node = node.getNext();
                k = k+1;
            end
            % after the step is made, increment the internal iterations counter of the object...
            obj.iter_ctr = obj.iter_ctr + 1;
            obj.prev_solution = res; % save solution y_n+1 as the previous solution
            res = gather(res);
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

close all; clear; clc; 

f = @(t) log(t+1) + sin(t) ;
f_t = @(t) 1/(t+1) + cos(t); 
coefs =  [23/12 ; -4/3; 5/12];
error = 0; 
tEnd = 1.5; 
dt = 1E-1;
attempt = 1; 
max_attempts = 6;

glob_err = zeros(max_attempts,3);

while attempt <= max_attempts
    clc; 
    display(['attempt Nr - ' num2str(attempt)]);
    dt = dt/2; 
    display(['dt = ' num2str(dt)]);
    nr_steps = 3; 

    %initialize the numerical solution with MS, the straightforward
    %imlpementation and the analytic solution for steps 0,1,2  

    %MSqueue with a full set of initial conditions provided!!! 
    init_conditions = [f_t(0),f_t(dt)];
    y_num_1 = [f(0),f(1*dt),f(2*dt)];
    solver_1 = MSqueueGPU(nr_steps,1,init_conditions,f(2*dt));

    %MSqueue with one initial condition provided and the rest self-computed
    y_num_2 = [f(0)];
    solver_2 = MSqueueGPU(nr_steps,1,[],f(0));

    y_num_2(2) = solver_2.make_step(f_t(0),dt); 
    y_num_2(3) = solver_2.make_step(f_t(dt),dt);

    y_num_out_of_the_box = [f(0),f(dt),f(2*dt)];
    y_analy = [f(0),f(dt),f(2*dt)];
    t = (nr_steps-1)*dt; 

    % y_num_solver(1) = solver.make_step(f_t(0),dt);
    % y_num_solver(2) = solver.make_step(f_t(dt),dt);
    % y_num_solver(3) = solver.make_step(f_t(2*dt),dt);
    ctr = 3; % starting with step 3...

    while t < tEnd

        ctr = ctr +1; 
        %solution at t+dt; = ctr*dt;
        y_num_1(ctr) = solver_1.make_step(f_t(t),dt);
        y_num_2(ctr) = solver_2.make_step(f_t(t),dt);
        y_num_out_of_the_box(ctr) = y_num_out_of_the_box(ctr-1) + dt*(f_t(t)*coefs(1) + f_t(t-dt)*coefs(2) +f_t(t-2*dt)*coefs(3));
        y_analy(ctr) = f(t+dt);

        if(mod(ctr,100) == 0)
            times = linspace(0,t,ctr);
            subplot(2,1,1);
            plot(times,y_num_1,'b',times,y_num_2,'g',times,y_num_out_of_the_box,'k',times,y_analy,'r');
            getframe;
        end

        t = t+dt; 
    end

    times = linspace(0,tEnd,ctr);
    error = zeros(ctr,3); 
    for k = 2:ctr
        error(k,1) = abs(y_num_1(k) - y_analy(k))/abs(y_analy(k));
        error(k,2) = abs(y_num_2(k) - y_analy(k))/abs(y_analy(k));
        error(k,3) = abs(y_num_out_of_the_box(k) - y_analy(k))/abs(y_analy(k));
    end
    
    glob_err(attempt,1) = max(error(:,1));
    glob_err(attempt,2) = max(error(:,2));
    glob_err(attempt,3) = max(error(:,3));
    attempt = attempt +1;
    
end
subplot(2,1,2); 
plot(1:max_attempts,glob_err(:,1),'b',1:max_attempts,glob_err(:,2),'g', 1:max_attempts,glob_err(:,3),'rx');

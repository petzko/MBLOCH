close all; clear; clc; 

f = @(t) log(t+1) + sin(t) ;
f_t = @(t) 1/(t+1) + cos(t); 
error = 0;  tEnd = 1.5; 
dt = 1E-1; attempt = 1; 
max_attempts = 8;
glob_err = zeros(max_attempts,3);

while attempt <= max_attempts
    
    clc; 
    display(['attempt Nr - ' num2str(attempt)]);
    dt = dt/2; 
    display(['dt = ' num2str(dt)]);
    nr_steps = 5; 

    %initialize the numerical solution with MS, the straightforward
    %imlpementation and the analytic solution for steps 0,1,2  

    %MS with a full set of initial conditions provided!!! 
    init_conditions = zeros(nr_steps-1,1);
    y_num_1 = [];
    for k = 1:nr_steps-1
        init_conditions(k) = f_t((k-1)*dt);
        y_num_1(k) = f((k-1)*dt);
    end
    solver_1 = MS_O1(nr_steps,1, init_conditions ,f((nr_steps-1)*dt));
    
    init_conditions = zeros(1,nr_steps-1);
    y_num_2 = [];
    for k = 1:nr_steps-1
        init_conditions(k) = f_t((k-1)*dt);
        y_num_2(k) = f((k-1)*dt);
    end
    solver_2 = MS(nr_steps,1, init_conditions ,f((nr_steps-1)*dt));
    
    
    
    
    y_analy = [f(0)];
    for k = 1:nr_steps-1
        y_analy(k+1) = f(k*dt);
    end
    t = (nr_steps-1)*dt;
    ctr = nr_steps-1; % starting with step 3...
    while t < tEnd

        ctr = ctr +1; 
        %solution at t+dt; = ctr*dt;
        y_num_1(ctr) = solver_1.make_step(f_t(t),dt);
        y_num_2(ctr) = solver_2.make_step(f_t(t),dt);

        y_analy(ctr) = f(t+dt);
        
        if(mod(ctr,100) == 0)
            times = linspace(0,t,ctr);
            plot(times,y_num_1,'b',times,y_num_2,'--g');
            getframe;
        end

        t = t+dt; 
    end

    times = linspace(0,tEnd,ctr);
    error = zeros(ctr,2); 
    for k = 2:ctr
        error(k,1) = abs(y_num_1(k) - y_analy(k));
        error(k,2) = abs(y_num_2(k) - y_analy(k));
    end
    
    glob_err(attempt,1) = max(error(:,1));
    glob_err(attempt,2) = max(error(:,2));
   
    attempt = attempt +1;
    
end
figure;
plot(1:max_attempts,glob_err(:,1),'b',1:max_attempts,glob_err(:,2),'g');

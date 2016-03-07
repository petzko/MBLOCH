close all; clear; clc; 
k = -1.0 +20*1i;

f = @(t) exp(k*t);
f_t = @(t) k*exp(k*t); 
coefs =  [23/12 ; -4/3; 5/12];
error = 0; 
tEnd = 10; 
dt = 1E-1;
attempt = 1; 
max_attempts = 8;

glob_err = zeros(max_attempts,3);

while attempt <= max_attempts
    
    display(['attempt Nr - ' num2str(attempt)]);
    dt = dt/2
    display(['dt = ' num2str(dt)]);
    nr_steps = 2; 

    %initialize the numerical solution with MS, the straightforward
    %imlpementation and the analytic solution for steps 0,1,2  

    %MSqueue with a full set of initial conditions provided!!! 
    init_conditions = [];
    y_num_1 = [];
    for k = 1:nr_steps-1
        init_conditions(k) = f_t((k-1)*dt);
        y_num_1(k) = f((k-1)*dt);
    end
    
    solver_1 = MS(nr_steps,1, init_conditions ,f((nr_steps-1)*dt));

    %MSqueue with one initial condition provided and the rest self-computed
    y_num_2 = [f(0)];
    y_analy = [f(0)];
    
    solver_2 = MS(nr_steps,1,[],f(0));
    for k = 1:nr_steps-1
        y_num_2(k+1) = solver_2.make_step(f_t((k-1)*dt),dt); 
        y_analy(k+1) = f(k*dt);
    end 
    
%     y_num_out_of_the_box = [f(0),f(dt),f(2*dt)];

    t = (nr_steps-1)*dt; 

    ctr = nr_steps-1; % starting with step 3...
    while t < tEnd

        ctr = ctr +1; 
        %solution at t+dt; = ctr*dt;
        y_num_1(ctr) = solver_1.make_step(f_t(t),dt);
        y_num_2(ctr) = solver_2.make_step(f_t(t),dt);

        y_analy(ctr) = f(t+dt);
        
% % %         if(mod(ctr,100) == 0)
% % %             times = linspace(0,t,ctr);
% % %             subplot(2,1,1);
% % %             plot(times,real(y_num_1),'b',times,real(y_num_2),'g',times,real(y_analy),'r');
% % % %             plot(times,y_num_1,'b',times,y_analy,'r');
% % % %             plot(times,y_num_2,'g',times,y_analy,'r');
% % %         
% % %                 getframe;
% % %         end

        t = t+dt; 
    end

    times = linspace(0,tEnd,ctr);
    error = zeros(ctr,3); 
    for k = 2:ctr
        error(k,1) = norm(y_num_1(k) - y_analy(k));
        error(k,2) = norm(y_num_2(k) - y_analy(k));
%         error(k,3) = abs(y_num_out_of_the_box(k) - y_analy(k))/abs(y_analy(k));
    end
    
    glob_err(attempt,1) = max(error(:,1));
    glob_err(attempt,2) = max(error(:,2));
    glob_err(attempt,3) = max(error(:,3));
    attempt = attempt +1;
    
end
figure;
subplot(2,1,2);
plot(1:max_attempts,glob_err(:,1),'b',1:max_attempts,glob_err(:,2),'g');

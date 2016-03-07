
%% initialize the solution

clear;clc;
f = @(t) log(t+1) + sin(t);
f_t = @(t) 1/(t+1) + cos(t); 
tEnd = 10; 
dt = 1E-2;
nr_steps = 3;
y_analy = 0;
ctr = 0;
t = 0;

%% setup the initial solution

while t <= tEnd
        ctr = ctr +1; 
        y_analy(ctr) = f(t);
        t = t+dt; 
end

%%  setup the CPU solver

init_conditions_cpu = [f_t(0),f_t(dt)];
y_num_cpu = zeros(length(y_analy),1);
y_num_cpu(1) = f(0);y_num_cpu(2) = f(1*dt);y_num_cpu(3) = f(2*dt);
solver_cpu = MSqueue(nr_steps,1,init_conditions_cpu,f(2*dt));
ctr_cpu = 3; % starting with step 3...
t_cpu = (nr_steps-1)*dt;
while t_cpu <= tEnd-dt

        ctr_cpu = ctr_cpu +1; 
        %solution at t+dt; = ctr*dt;
        y_num_cpu(ctr_cpu) = solver_cpu.make_step(f_t(t_cpu),dt);
        t_cpu = t_cpu+dt; 
end


%% setup the GPU solver

init_conditions_gpu = [f_t(0),f_t(dt)];
y_num_gpu = zeros(length(y_analy),1);
y_num_gpu(1)  = f(0); y_num_gpu(2) = f(1*dt); y_num_gpu(3) = f(2*dt);
y_num_gpu = gpuArray(y_num_gpu);
solver_gpu = MSqueueGPU(nr_steps,1,init_conditions_gpu,f(2*dt));
ctr_gpu = 3; % starting with step 3...
t_gpu = (nr_steps-1)*dt;
while t_gpu < tEnd-dt
        ctr_gpu = ctr_gpu +1; 
        %solution at t+dt; = ctr*dt;
        y_num_gpu(ctr_gpu) = solver_gpu.make_step(f_t(t_gpu),dt);
        t_gpu = t_gpu+dt; 
end
plot(y_num_gpu)

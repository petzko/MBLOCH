
%initial condition
f_0 = 100; 
%analytical solution
f = @(t_) f_0*exp(-t_/15);

% differential equaiton to solve
df_dt = @(t_) -1/15*f(t_);

%number of time steps to propagate
N_t = 1000;
%duration of each single time-step
dt = 1e-2; 

% numerical & analytical solution vectors
f_num = zeros(N_t,1);
f_analytical = zeros(N_t,1);

%for Olfa -> set nr_step = 5
nr_steps = 5;
solver = MS(nr_steps,N,[],f_0)


t = dt; 
ctr = 1; 
while(t < N_t*dt)
    % compute the analytical solution
    f_analytical(ctr) =  f(t);
    f_num(ctr) = solver.make_step(df_dt(t-dt),dt);
    
     if mod(ctr,1) == 0
        plot(1:N_t,f_analytical,'-r',1:N_t,f_num,'--b')
        getframe;
    end
    
    ctr = ctr + 1; 
    t = t+dt;
end





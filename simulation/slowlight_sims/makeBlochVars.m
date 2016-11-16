function dat = makeBlochVars(settings,dat)

dtype = dat.dtype;


%%%%  gain section
%populations
dat.r22 = 0.9*ones(settings.N,1,dtype);
dat.r11 = 1-dat.r22;
%coherence
% dat.n21 =  1E-15*((rand(settings.N,1,dtype)-0.5) +1i*(rand(settings.N,1,dtype)-0.5));
dat.n21 =  1E-20*((ones(settings.N,1,dtype)-0.5) +1i*(ones(settings.N,1,dtype)-0.5));
dat.n21_t = 0*dat.n21;
%solvers
dat.r22_solver = MS(settings.nr_steps,settings.N,[],dat.r22);
dat.r11_solver = MS(settings.nr_steps,settings.N,[],dat.r11);
dat.n21_solver = MS(settings.nr_steps,settings.N,[],dat.n21);

%%%% slow light section
%populations
dat.r_ee = 0*ones(settings.N,1,dtype);
dat.r_gg = ones(settings.N,1,dtype); 
dat.r_ss = 0*ones(settings.N,1,dtype);

% coherence terms
dat.r_se = 1E-20*((rand(settings.N,1,dtype)-0.5) +1i*(rand(settings.N,1,dtype)-0.5)); 
dat.n_eg = 1E-20*((rand(settings.N,1,dtype)-0.5) +1i*(rand(settings.N,1,dtype)-0.5));
dat.n_sg = 1E-20*((rand(settings.N,1,dtype)-0.5) +1i*(rand(settings.N,1,dtype)-0.5));

dat.r_se = 1E-20*((ones(settings.N,1,dtype)-0.5) +1i*(ones(settings.N,1,dtype)-0.5)); 
dat.n_eg = 1E-20*((ones(settings.N,1,dtype)-0.5) +1i*(ones(settings.N,1,dtype)-0.5));
dat.n_sg = 1E-20*((ones(settings.N,1,dtype)-0.5) +1i*(ones(settings.N,1,dtype)-0.5));

dat.r_ss_solver = MS(settings.nr_steps,settings.N,[],dat.r_ss); 
dat.r_gg_solver = MS(settings.nr_steps,settings.N,[],dat.r_gg);
dat.r_ee_solver = MS(settings.nr_steps,settings.N,[],dat.r_ee); 

dat.n_sg_solver = MS(settings.nr_steps,settings.N,[],dat.n_sg);
dat.r_se_solver = MS(settings.nr_steps,settings.N,[],dat.r_se); 
dat.n_eg_solver = MS(settings.nr_steps,settings.N,[],dat.n_eg);

dat.P = zeros(settings.N,1); 
dat.P_t = zeros(settings.N,1); 


end
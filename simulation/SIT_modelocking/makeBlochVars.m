function dat = makeBlochVars(settings,dat)
    
    dtype = dat.dtype;
    
    
    %%%%  gain section 
    %populations
    dat.r22g = ones(dat.Ng,1,dtype); 
    dat.r11g = zeros(dat.Ng,1,dtype);
    %coherence 
    dat.n21g = 0*1E-15*((rand(dat.Ng,1,dtype)-0.5) +1i*(rand(dat.Ng,1,dtype)-0.5)); 
    dat.n21g_t = 0*dat.n21g;
    %solvers
    dat.r22g_solver = MS(settings.nr_steps,dat.Ng,[],dat.r22g); 
    dat.r11g_solver = MS(settings.nr_steps,dat.Ng,[],dat.r11g);
    dat.n21g_solver = MS(settings.nr_steps,dat.Ng,[],dat.n21g); 
     
    %%%% absorber section 
    %populations

%     r22a -> ull
     dat.r22a = zeros(dat.Na,1,dtype); 
     dat.r11a = ones(dat.Na,1,dtype);
    

    %cohrence
    dat.n21a =  0*1E-15*((rand(dat.Na,1,dtype)-0.5)*0 +1i*(rand(dat.Na,1,dtype)-0.5)); 
    dat.n21a_t = 0*dat.n21a;
    
    %solvers
    dat.r22a_solver = MS(settings.nr_steps,dat.Na,[],dat.r22a); 
    dat.r11a_solver = MS(settings.nr_steps,dat.Na,[],dat.r11a);
    dat.n21a_solver = MS(settings.nr_steps,dat.Na,[],dat.n21a); 
   
    %global coherence and dervative vars
    dat.n21 = zeros(settings.N,1,dat.dtype); 
    dat.n21_t = zeros(settings.N,1,dat.dtype); 
    
    
end
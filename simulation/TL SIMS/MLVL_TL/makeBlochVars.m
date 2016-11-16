function dat = makeBlochVars(settings,dat)
    
    dtype = dat.dtype;
    

    tau2B = 1E-12/dat.W(dat.LLL,dat.DEPOP);
    tau3B = 1E-12/dat.W(dat.ULL,dat.DEPOP); 
    Gamma0 = dat.i0*1E3/(settings.Ltot*1E-3*settings.Lp*1E-9*dat.Ncarriers*Constants('q0'));
    init_r33 = tau2B*tau3B/(tau2B-tau3B)*(Gamma0-1/tau2B);
        
    init_r22 = 1-init_r33;
    % population vectors together with first derivative vectors
    dat.r110 = dat.trace_rho*zeros(settings.N,1,dtype); 
    dat.r330 = dat.trace_rho*init_r33*ones(settings.N,1,dtype);
    dat.r220 = dat.trace_rho*init_r22*ones(settings.N,1,dtype);
    dat.populations = zeros(settings.N,dat.N_rest,dtype);

    dat.rate_eqn_solvers = [];
    for p=1:dat.N_rest
        dat.rate_eqn_solvers{p} =  MS(settings.nr_steps,settings.N,[],dat.populations(:,p));
    end
 
    % coherence terms together with first derivative vectors!
    dat.r130 =0*((ones(settings.N,1,dtype)-0.5) +1i*(ones(settings.N,1,dtype)-0.5)); 
   
    
    dat.n32p =  1E-15*((rand(settings.N,1,dtype)-0.5) +1i*(rand(settings.N,1,dtype)-0.5));
    dat.n32m =  1E-15*((rand(settings.N,1,dtype)-0.5) +1i*(rand(settings.N,1,dtype)-0.5)); 
    
    dat.n12p =  1E-15*((rand(settings.N,1,dtype)-0.5) +1i*(rand(settings.N,1,dtype)-0.5)); 
    dat.n12m =  1E-15*((rand(settings.N,1,dtype)-0.5) +1i*(rand(settings.N,1,dtype)-0.5)); 
 
    
    dat.r110_solver = MS(settings.nr_steps,settings.N,[],dat.r110); 
    dat.r330_solver = MS(settings.nr_steps,settings.N,[],dat.r330);
    dat.r220_solver = MS(settings.nr_steps,settings.N,[],dat.r220);
    
    %the solvers for the rest of the auxilliary rates.
    dat.n12p_solver = MS(settings.nr_steps,settings.N,[],dat.n12p); 
    dat.n12m_solver = MS(settings.nr_steps,settings.N,[],dat.n12m);
    dat.n32p_solver = MS(settings.nr_steps,settings.N,[],dat.n32p);
    dat.n32m_solver = MS(settings.nr_steps,settings.N,[],dat.n32m);
    dat.r130_solver = MS(settings.nr_steps,settings.N,[],dat.r130);
    
   
    
    if (settings.shb > 0)
        
        dat.r11p = zeros(settings.N,1,dtype);
        dat.r33p = zeros(settings.N,1,dtype);
        dat.r22p = zeros(settings.N,1,dtype);
        dat.r13p = zeros(settings.N,1,dtype); 
        dat.r13m = zeros(settings.N,1,dtype);
        
        dat.r13p_solver = MS(settings.nr_steps,settings.N,[],dat.r13p); 
        dat.r13m_solver = MS(settings.nr_steps,settings.N,[],dat.r13m);
        
        dat.r11p_solver =MS(settings.nr_steps,settings.N,[],dat.r11p);
        dat.r33p_solver = MS(settings.nr_steps,settings.N,[],dat.r33p);
        dat.r22p_solver = MS(settings.nr_steps,settings.N,[],dat.r22p);
    
    end
 

end
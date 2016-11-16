function dat = makeBlochVars(settings,dat)
       
    dtype = dat.dtype;
    

    % population vectors together with first derivative vectors
    dat.r110 = 1/3*ones(settings.N,1,dtype); 
    dat.r330 = 1/3*ones(settings.N,1,dtype);
    dat.r220 = 1/3*ones(settings.N,1,dtype);
    dat.rRES = zeros(settings.N,1,dtype); 
    
    % coherence terms together with first derivative vectors!
    dat.r130 =1E-15*((rand(settings.N,1,dtype)) +1i*(rand(settings.N,1,dtype)));
   
    dat.n32p =  1E-15*((rand(settings.N,1,dtype)) +1i*(rand(settings.N,1,dtype)));
    dat.n32m = 1E-15*((rand(settings.N,1,dtype)) +1i*(rand(settings.N,1,dtype)));
    
    dat.n12p =  1E-15*((rand(settings.N,1,dtype)) +1i*(rand(settings.N,1,dtype))); 
    dat.n12m =  1E-15*((rand(settings.N,1,dtype)) +1i*(rand(settings.N,1,dtype)));; 
 
    
    dat.r110_solver = MS(settings.nr_steps,settings.N,[],dat.r110); 
    dat.r330_solver = MS(settings.nr_steps,settings.N,[],dat.r330);
    dat.r220_solver = MS(settings.nr_steps,settings.N,[],dat.r220);
    dat.rRES_solver = MS(settings.nr_steps,settings.N,[],dat.rRES);

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
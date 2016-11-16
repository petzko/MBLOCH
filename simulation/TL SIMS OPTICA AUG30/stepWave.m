function dat = stepWave(settings,dat)

    %%%%%%%%%%%%%%%%%%%%%%%%
    %interaction strength factor -trace_rho times dipR;
    factor = dat.factor*dat.dipR;
    dat.U = dat.U_solver.make_step(factor.*dat.n32p,factor.*dat.n32p_t,dat.losses,dat.dt);
    dat.V = dat.V_solver.make_step(factor.*dat.n32m,factor.*dat.n32m_t,dat.losses,dat.dt);
    
    %set the boundaries... and obtain the final solution
    dat.U = dat.U_solver.set_bdry(dat.V(1),'no');
    dat.V = dat.V_solver.set_bdry('no',dat.U(end));
    %%%%%%%%%%%%%%%%%%%%%%%%


end
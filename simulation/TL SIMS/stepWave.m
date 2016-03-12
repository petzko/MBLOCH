function dat = stepWave(settings,dat)

    %%%%%%%%%%%%%%%%%%%%%%%%
    
    dat.U = dat.U_solver.make_step(-1i*dat.c*dat.n32p,-1i*dat.c*dat.n32p_t,dat.losses,dat.dt);
    dat.V = dat.V_solver.make_step(-1i*dat.c*dat.n32m,-1i*dat.c*dat.n32m_t,dat.losses,dat.dt);
    
    %set the boundaries... and obtain the final solution
    dat.U = dat.U_solver.set_bdry(dat.V(1),'no');
    dat.V = dat.V_solver.set_bdry('no',dat.U(end));
    %%%%%%%%%%%%%%%%%%%%%%%%


end
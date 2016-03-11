function dat = stepWave(settings,dat)

    %%%%%%%%%%%%%%%%%%%%%%%%
    losses = -dat.c*dat.l_0.*ones(settings.N,1);
    
    dat.U = dat.U_solver.make_step(-1i*dat.c*dat.n32p,-1i*dat.c*dat.n32p_t,losses,dat.dt);
    dat.V = dat.V_solver.make_step(-1i*dat.c*dat.n32m,-1i*dat.c*dat.n32m_t,losses,dat.dt);
    
    %set the boundaries... and obtain the final solution
    dat.U = dat.U_solver.set_bdry(dat.V(1),'no');
    dat.V = dat.V_solver.set_bdry('no',dat.U(end));
    %%%%%%%%%%%%%%%%%%%%%%%%


end
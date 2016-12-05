function dat = stepWave(dat,P,P_t,M,M_t,losses)

    %%%%%%%%%%%%%%%%%%%%%%%%
    %interaction strength factor -trace_rho times dipR;
    factor =-1i*dat.c;
    dat.U = dat.U_solver.make_step(factor.*P,factor.*P_t,-dat.c*losses,dat.dt);
    dat.V = dat.V_solver.make_step(factor.*M,factor.*M_t,-dat.c*losses,dat.dt);
    
    %set the boundaries... and obtain the final solution
    dat.U = dat.U_solver.set_bdry(dat.V(1),'no');
    dat.V = dat.V_solver.set_bdry('no',dat.U(end));
    %%%%%%%%%%%%%%%%%%%%%%%%


end
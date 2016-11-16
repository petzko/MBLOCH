function dat = stepWave(settings,dat)

    %interaction strength factor;
    ftr = dat.factor.*dat.dipR.*dat.i_core;
    dat.U = dat.wave_solver.make_step(ftr.*dat.P,ftr.*dat.P_t,-dat.c*dat.losses,dat.dt);
    %set the boundaries... and obtain the final solution
    dat.U = dat.wave_solver.set_bdry(dat.message(dat.t),'no');
    %%%%%%%%%%%%%%%%%%%%%%%%

end
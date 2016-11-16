function dat = stepWave(settings,dat)

    %interaction strength factor;
    ftr = dat.factor;
    
    %total polarization -> combination of gain and abosrption polarizations
    dat.n21(dat.gain_start:dat.gain_end) = dat.n21g;
    dat.n21_t(dat.gain_start:dat.gain_end) = dat.n21g_t;
    dat.n21(dat.abs_start:dat.abs_end) = dat.n21a;
    dat.n21_t(dat.abs_start:dat.abs_end) = dat.n21a_t;

    dat.U = dat.wave_solver.make_step(dat.core.*ftr.*dat.n21,...
        dat.core.*ftr.*dat.n21_t,dat.core.*dat.losses,dat.dt);
    %set the boundaries... and obtain the final solution
    tp = dat.tp;
    A0 = dat.A0;
    t0 = dat.t0;
    
    dat.U = dat.wave_solver.set_bdry(A0*sech((dat.t-t0)/tp)+ ...
        dat.U(end)*settings.R,'no');
    %%%%%%%%%%%%%%%%%%%%%%%%

end
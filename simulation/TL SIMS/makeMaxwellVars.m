function dat = makeMaxwellVars(settings,dat)

    x_0 = settings.Ltot/7; tp =1;
    aE_in = @(z,time) exp(-(time-(z-x_0)/dat.c).^2/tp^2);
    
    dat.U = aE_in(dat.x,0);
    %normalize
    ampl = 3;
    dat.U = ampl*dat.U/max(abs(dat.U));
    dat.V = 0*dat.U;
    
    dat.U_solver = RNFDSolver(settings.N,dat.dx,+1,dat.c, dat.U);
    dat.V_solver = RNFDSolver(settings.N,dat.dx,-1,dat.c,dat.V);
    

end
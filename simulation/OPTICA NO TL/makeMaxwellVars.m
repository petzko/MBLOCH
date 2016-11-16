function dat = makeMaxwellVars(settings,dat)

    x_0 = 0*settings.Ltot/7; tp =1;
    aE_in = @(z,time) exp(-(time-(z-x_0)/dat.c).^2/tp^2);
    
    dat.U = aE_in(dat.x,0);
    if (strcmp(dat.dtype,'single'))
        dat.U = single(dat.U);
    end
    %normalize
    ampl = 1;
    dat.U = ampl*dat.U/max(abs(dat.U));
    dat.V = 0*dat.U;
 
    dat.U_solver = RNFDSolver(settings.N,dat.dx,+1,dat.c, dat.U);
    dat.V_solver = RNFDSolver(settings.N,dat.dx,-1,dat.c,dat.V);
    
    dat.losses = -dat.c*dat.l_0.*ones(settings.N,1,dat.dtype);
    dat.factor = -1i*dat.c*dat.trace_rho;
    
end
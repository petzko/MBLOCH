function dat = makeMaxwellVars(settings,dat)

      
    dat.U = zeros(settings.N,1);
    if (strcmp(dat.dtype,'single'))
        dat.U = single(dat.U);
    end
    % approximate magnitude
    magn = 1e-5;
    %electric field vector
    dat.U = magn*( (rand(settings.N,1,dat.dtype)-.5)+ ... 
                   1i*(rand(settings.N,1,dat.dtype) - .5));
    %electric field solver 
    dat.wave_solver = RNFDSolver(settings.N,dat.dx,+1,dat.c, dat.U);
    
    g1 = dat.gain_start; g2 = dat.gain_end; % g2-g1 should =  dat.Ng-1 
    a1 = dat.abs_start; a2 = dat.abs_end; % a2-a1 should = dat.Na-1
    
    %position dependent linear losses both in the gain and absorption
    %sections 
    dat.losses = ones(settings.N,1,dat.dtype);
    dat.losses(g1:g2) = -dat.c*dat.lg
    dat.losses(a1:a2) = -dat.c*dat.la;

    %position dependent normalization constant (trace)!!!
    dat.factor = ones(settings.N,1,dat.dtype); 
    dat.factor(g1:g2) = -1i*dat.c*dat.trace_rho_g;
    dat.factor(a1:a2) = -1i*dat.c*dat.trace_rho_a;
    dat.factor = dat.factor.*dat.dipR;
end
function dat = makeMaxwellVars(dat)
    dat.U = zeros(dat.N,1); 
    dat.V = dat.U;
    
    dat.U_solver = RNFDSolver(dat.N,dat.dx,+1,dat.c, dat.U);
    dat.V_solver = RNFDSolver(dat.N,dat.dx,-1,dat.c, dat.V);
    
    
end
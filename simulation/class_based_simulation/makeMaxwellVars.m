function dat = makeMaxwellVars(params,dat)
    dat.U = zeros(params.N,1); 
    dat.V = dat.U;
    
    dat.U_solver = RNFDSolver(params.N,params.dx,+1,params.c, dat.U);
    dat.V_solver = RNFDSolver(params.N,params.dx,-1,params.c,dat.V);
    
    
end
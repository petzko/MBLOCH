function dat = updateBloch(settings,dat)

    %%%%%%%%%%%%%%%%%%%%%%%%
   
    
    dat.r110 = dat.r110_solver.get_latest_solution();
    dat.r330 = dat.r330_solver.get_latest_solution();
    dat.r220 = dat.r220_solver.get_latest_solution();

    
    dat.r130 = dat.r130_solver.get_latest_solution(); dat.n32p = dat.n32p_solver.get_latest_solution();
    dat.n32m = dat.n32m_solver.get_latest_solution(); dat.n12p = dat.n12p_solver.get_latest_solution(); 
    dat.n12m = dat.n12m_solver.get_latest_solution();
    
    if(settings.shb > 0)
        dat.r11p = dat.r11p_solver.get_latest_solution();
        dat.r33p = dat.r33p_solver.get_latest_solution();
        dat.r22p = dat.r22p_solver.get_latest_solution();
        dat.r13p = dat.r13p_solver.get_latest_solution();
        dat.r13m = dat.r13m_solver.get_latest_solution();
    end


end
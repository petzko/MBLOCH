function dat = updateBloch(settings,dat)

    %%%%%%%%%%%%%%%%%%%%%%%%
    dat.r22g = dat.r22g_solver.get_latest_solution();
    dat.r11g = dat.r11g_solver.get_latest_solution();
    dat.r22a = dat.r22a_solver.get_latest_solution();
    dat.r11a = dat.r11a_solver.get_latest_solution();
   
    dat.n21g = dat.n21g_solver.get_latest_solution();
    dat.n21a = dat.n21a_solver.get_latest_solution();
   

end
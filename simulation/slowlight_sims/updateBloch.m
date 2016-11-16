function dat = updateBloch(settings,dat)

dat.r_ss = dat.r_ss_solver.get_latest_solution();
dat.r_ee = dat.r_ee_solver.get_latest_solution();
dat.r_gg = dat.r_gg_solver.get_latest_solution();


dat.r_se = dat.r_se_solver.get_latest_solution();
dat.n_eg = dat.n_eg_solver.get_latest_solution();
dat.n_sg = dat.n_sg_solver.get_latest_solution();


%%%%%%%%%%%%%%%%%%%%%%%%
dat.r22 = dat.r22_solver.get_latest_solution();
dat.r11 = dat.r11_solver.get_latest_solution();
dat.n21 = dat.n21_solver.get_latest_solution();


end
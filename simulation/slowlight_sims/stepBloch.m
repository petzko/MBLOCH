function dat = stepBloch(settings,dat)


%%%%% IMPLEMENTATION 2
%%%% GAIN SECTION
%populations
lmInteraction = conj(dat.U).*dat.n21;
r22_t = 1i/2*(lmInteraction-conj(lmInteraction)) -...
    (dat.r22-dat.r22_0)/dat.T1;
dat.r22_solver.make_step(r22_t,dat.dt);

r11_t = -1i/2*(lmInteraction-conj(lmInteraction)) - ...
    (dat.r11-dat.r11_0)/dat.T1;
dat.r11_solver.make_step(r11_t,dat.dt);

%coherences
n21_t = dat.dE21.*dat.n21 + ...
    1i/2*dat.dipR.*dat.U.*(dat.r22-dat.r11);
dat.n21_solver.make_step(dat.n21_t,dat.dt);


%%%%  SLOW DOWN SECTION
W = dat.W; Ex_ = dat.Ex_; Spin_ = dat.Spin_; Grnd_ = dat.Grnd_; G = dat.G;
i_dipR = dat.dipR_inv; %inverse dipole ratio
%%%% POPULATIONS
rss_t = 1i*dat.O_se*(dat.r_se-conj(dat.r_se)) +( W(Ex_,Spin_)  )*dat.r_ee +...
    W(Grnd_,Spin_)*dat.r_gg - G(Spin_)*dat.r_ss ;
dat.r_ss_solver.make_step(rss_t,dat.dt);

lmInteraction = i_dipR.*conj(dat.U).*dat.n_eg;
ree_t = -1i*dat.O_se*(dat.r_se-conj(dat.r_se)) + ...
    1i/2.*(lmInteraction-conj(lmInteraction)) + dat.r_ss*(W(Spin_,Ex_)) +...
    dat.r_gg*W(Grnd_,Ex_) - G(Ex_)*dat.r_ee;
dat.r_ee_solver.make_step(ree_t,dat.dt);

rgg_t = -1i/2.*(lmInteraction-conj(lmInteraction)) + ...
    dat.r_ss*(W(Spin_,Grnd_)) + dat.r_ee*W(Ex_,Grnd_) - G(Grnd_)*dat.r_gg;
dat.r_gg_solver.make_step(rgg_t,dat.dt);

%%% coherences! r13 n32 n12
rse_t = dat.dE_se*dat.r_se + 1i*dat.O_se*(dat.r_ss - dat.r_ee) + ...
    1i/2*i_dipR.*conj(dat.U).*dat.n_sg;
dat.r_se_solver.make_step(rse_t,dat.dt);

neg_t = dat.dE_eg*dat.n_eg + 1i/2*i_dipR.*dat.U.*(dat.r_ee-dat.r_gg) - ...
    1i*dat.O_se*dat.n_sg;
dat.n_eg_solver.make_step(neg_t,dat.dt);

n12_t = dat.dE_sg*dat.n_sg + 1i/2*i_dipR.*dat.U.*dat.r_se - 1i*dat.O_se*dat.n_eg;
dat.n_sg_solver.make_step(n12_t,dat.dt);

dat.P = dat.n_eg.*dat.core_SD+dat.n21.*dat.core_GAIN;
dat.P_t= neg_t.*dat.core_SD+n21_t.*dat.core_GAIN;

end
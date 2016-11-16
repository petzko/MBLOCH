function dat = stepBloch(settings,dat)
    
    % integer indices specifying the points on the grid where the 
    % gain and absorption secitons start and end respectively 
    g1 = dat.gain_start; g2 = dat.gain_end; % g2-g1 should =  dat.Ng-1 
    a1 = dat.abs_start; a2 = dat.abs_end; % a2-a1 should = dat.Na-1
    % both the dipole moment ratio and the detuning WILL be position
    % dependent, dependign on which section (gain or absorption) they lie
    % withi. 
  

%%%%% IMPLEMENTATION  
    %%%% GAIN SECTION
    %populations
    lmInteraction = conj(dat.U(g1:g2)).*dat.n21g;
    r22g_t = 1i/2*(lmInteraction-conj(lmInteraction)) -...
                    (dat.r22g-settings.r22g_0)/settings.T1_g;
                
    dat.r22g_solver.make_step(r22g_t,dat.dt);
    r11g_t = -1i/2*(lmInteraction-conj(lmInteraction)) - ...
                    (dat.r11g-settings.r11g_0)/settings.T1_g;
    dat.r11g_solver.make_step(r11g_t,dat.dt);
    
    %coherences
    dat.n21g_t = dat.dE21(g1:g2).*dat.n21g + ...
                    1i/2*dat.dipR(g1:g2).*dat.U(g1:g2).*(dat.r22g-dat.r11g);
    dat.n21g_solver.make_step(dat.n21g_t,dat.dt);
   
    
    %%%%  ABSORPTION SECTION
    %populations
    lmInteraction = dat.dipR_inv(a1:a2).*conj(dat.U(a1:a2)).*dat.n21a;
    r22a_t = 1i/2*(lmInteraction-conj(lmInteraction)) - ...
                    (dat.r22a-settings.r22a_0)/settings.T1_a;
    dat.r22a_solver.make_step(r22a_t,dat.dt);
    r11a_t = -1i/2*(lmInteraction-conj(lmInteraction)) - ...
                    (dat.r11a-settings.r11a_0)/settings.T1_a;
    dat.r11a_solver.make_step(r11a_t,dat.dt);
    
    %coherences
    dat.n21a_t = dat.dE21(a1:a2).*dat.n21a + ...
                    1i/2*dat.dipR_inv(a1:a2).*dat.U(a1:a2).*(dat.r22a-dat.r11a);
    dat.n21a_solver.make_step(dat.n21a_t,dat.dt);

   
end
function dat = stepTransLine(settings,dat)
%%%%% Begin TRANSMISSION LINE EQUATIONS %%%%%%%%%%%%%

dat.J_TL = (Constants('q0')*(settings.Lp*1E-9)*dat.Ncarriers*1E12*dat.rates)/1E6; %in A/mm^2


%% Transmission Line Voltage wave equation
% dat.J_TL_t = (Constants('q0')*(settings.Lp*1E-9)*dat.Ncarriers*1E12*dat.rates_t)/1E6; %in A/mm^2
% dat.v_tmp(2:end-1) = 2*dat.v_TL(2:end-1)-dat.v_TL_old(2:end-1) + dat.Acoeff*(dat.v_TL(3:end)-2*dat.v_TL(2:end-1)+dat.v_TL(1:end-2)) - dat.Bcoeff*dat.J_TL_t(2:end-1,1);
% dat.v_tmp(1) = dat.v0_t(dat.t);
%%%%% END Transmission Line Voltage wave equation %%%%%%%%%%%%%


%% Transmission Line Current-Voltage equation
%     dat.i_TL(2:end) = dat.i_TL(2:end)*(1-dat.dt*dat.R_Au)-dat.Dcoeff*(dat.v_TL(2:end)-dat.v_TL(1:end-1));
%     dat.i_TL(1) = dat.i0_t(dat.t);
%     dat.v_tmp(1:end-1) = dat.v_TL(1:end-1)-dat.Ecoeff*(dat.i_TL(2:end)-dat.i_TL(1:end-1))-dat.Fcoeff*dat.J_TL(1:end-1);
%%%%% END Transmission Line Current-Voltage equation %%%%%%%%%%%%%

%     ABSBDRY= dat.v_TL(end-1)+(settings.nTHz-settings.nRF)/(settings.nTHz+settings.nRF)*(dat.v_tmp(end-1)-dat.v_TL(end));
%     OCBDRY =  2*dat.v_TL(end)-dat.v_TL_old(end)+2*dat.Acoeff*(dat.v_TL(end-1)-dat.v_TL(end))-0*dat.Bcoeff*dat.J_TL_t(end,1);
%     dat.v_tmp(end) = OCBDRY;

%i->v (unstable)
%     dat.v_tmp(2:end) = dat.v_TL(2:end)-dat.Ecoeff*(dat.i_TL(2:end)-dat.i_TL(1:end-1))-dat.Fcoeff*dat.J_TL(2:end);
%     dat.v_tmp(1) =(dat.u_t(dat.t)-dat.Z_kV_A*dat.i_TL(1)*dat.width_mm )/dat.Lp_mm;
%     
%     dat.i_TL(1:end-1) = dat.i_TL(1:end-1)-dat.Dcoeff*(dat.v_TL(2:end)-dat.v_TL(1:end-1));
%     dat.i_TL(end) = 0;

%v->i (stable)
  
    %first solve and set the voltage! 
    dat.v_tmp(1:end-1) = dat.v_TL(1:end-1)-dat.Ecoeff*(dat.i_TL(2:end)-dat.i_TL(1:end-1))-dat.Fcoeff*dat.J_TL(1:end-1);
    OCBDRY =  2*dat.v_TL(end)-dat.v_TL_old(end)+2*dat.Acoeff*(dat.v_TL(end-1)-dat.v_TL(end));
    dat.v_tmp(end) = OCBDRY;
    dat.v_TL_old = dat.v_TL; dat.v_TL =  dat.v_tmp ;
    
    %then solve and set the current
    dat.i_TL(2:end) = dat.i_TL(2:end)-dat.Dcoeff*(dat.v_TL(2:end)-dat.v_TL(1:end-1));
    dat.i_TL(1) = (dat.u_t(dat.t)-dat.v_TL(1)*dat.Lp_mm)/(dat.Z_kV_A*dat.width_mm);


%     display(num2str([dat.v_TL(1), dat.i_TL(1)]))

end
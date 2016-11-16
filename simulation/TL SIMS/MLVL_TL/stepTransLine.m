function dat = stepTransLine(settings,dat)
 %%%%% Begin TRANSMISSION LINE EQUATIONS %%%%%%%%%%%%%
   
    W = dat.W; INJ = dat.INJ; ULL =dat.ULL; LLL = dat.LLL; DEPOP =dat.DEPOP;
    dat.J_TL = (Constants('q0')*(settings.Lp*1E-9)*dat.Ncarriers/dat.trace_rho*1E12*dat.rates)/1E6; %in A/mm^2

    
    %% Transmission Line Voltage wave equation
%     dat.J_TL_t = (Constants('q0')*(settings.Lp*1E-9)*dat.Ncarriers/dat.trace_rho*1E12*dat.rates_t)/1E6; %in A/mm^2
%     dat.v_tmp(2:end-1) = 2*dat.v_TL(2:end-1)-dat.v_TL_old(2:end-1) + dat.Acoeff*(dat.v_TL(3:end)-2*dat.v_TL(2:end-1)+dat.v_TL(1:end-2)) - dat.Bcoeff*dat.J_TL_t(2:end-1); 
%     dat.v_tmp(1) = dat.v0_t(dat.t); 
%     %open circuit bdry cond. (reflection)
%     %explicit:  dat.v_tmp(end) = 2*dat.v_TL(end)-dat.v_TL_old(end)+2*dat.Acoeff*(dat.v_TL(end-1)-dat.v_TL(end))- dat.Bcoeff*dat.J_TL_t(end); 
%     %implicit:  dat.v_tmp(end) = 1/(1+2*Acoeff)*(2*v_TL_new(end)-v_TL_old(end)+2*Acoeff*v_tmp(end-1)); 
%     %absorbing bdry cond:
%     dat.v_tmp(end)= dat.v_TL(end-1)+(settings.nTHz-settings.nRF)/(settings.nTHz+settings.nRF)*(dat.v_tmp(end-1)-dat.v_TL(end));
%     dat.v_TL_old = dat.v_TL; dat.v_TL = v_tmp ;  
%     
    %%%%% END Transmission Line Voltage wave equation %%%%%%%%%%%%%

    %% Transmission Line Current-Voltage equation
    
    dat.i_TL(2:end) = dat.i_TL(2:end)*(1-dat.dt*dat.R_Au)-dat.Dcoeff*(dat.v_TL(2:end)-dat.v_TL(1:end-1)); 
    dat.i_TL(1) = dat.i0_t(dat.t); 
    dat.v_tmp(1:end-1) = dat.v_TL(1:end-1)-dat.Ecoeff*(dat.i_TL(2:end)-dat.i_TL(1:end-1))-dat.Fcoeff*dat.J_TL(1:end-1);
    
    
    %set bdry conditions for v at right end! 
      %open circuit bdry cond. (reflection)
%     %explicit:  
%      dat.v_tmp(end) = 2*dat.v_TL(end)-dat.v_TL_old(end)+2*(settings.nTHz/settings.nRF)^2*(dat.v_TL(end-1)-dat.v_TL(end)); 
     dat.v_tmp(end) =  dat.v_tmp(end-1);
    %absorbing boundaries. 
%     %absorbing bdry cond:
%     dat.v_tmp(end)= dat.v_TL(end-1)+(settings.nTHz-settings.nRF)/(settings.nTHz+settings.nRF)*(dat.v_tmp(end-1)-dat.v_TL(end));

    %%%%% END Transmission Line Current-Voltage equation %%%%%%%%%%%%%
     dat.v_TL_old = dat.v_TL; dat.v_TL =  dat.v_tmp ;  
   
end
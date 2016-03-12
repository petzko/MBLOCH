function dat = makeTransLineVars(settings,dat)

    %make transline vars

    INJ = dat.INJ; ULL = dat.ULL; LLL = dat.LLL; DEPOP = dat.DEPOP; G = dat.G;
    %%% curr density calculation:
    if settings.Zorder > 0
        dat.rates = dat.r110.*dat.W(INJ,DEPOP) + dat.r220.*dat.W(LLL,DEPOP) + (dat.r330).*dat.W(ULL,DEPOP);
    else
        dat.rates = dat.r110.*dat.W_fit{INJ,DEPOP}(v_TL) + (dat.r220).*dat.W_fit{LLL,DEPOP}(v_TL) + (dat.r330).*dat.W_fit{ULL,DEPOP}(v_TL);
    end

    %convert the curr density in A/m^2 into A/mm^2 by dividing by 1E6;
    dat.J_TL = (Constants('q0')*(settings.Lp*1E-9)*dat.Ncarriers/dat.trace_rho*1E12*dat.rates)/1E6; %in A/mm^2
    dat.J_TL_t = 0*dat.J_TL;
    eps_ch = 1E15*Constants('eps0'); mu_ch  = 1E3*Constants('mu0');

    %impedance:
    dat.Z0 = sqrt(mu_ch/eps_ch)/settings.nRF;

    %voltage
    dat.v0 = settings.bias/1E1;
    %voltage modulation funciton
    dat.v0_t = @(tm) dat.v0*(1+settings.modA*sin(2*pi*settings.modF*dat.f_R*tm));
    dat.v_TL = dat.v0*ones(settings.N,1,dat.dtype); % transmission line voltage per unit length (in units kV/mm);
    dat.v_TL_old = dat.v0*ones(settings.N,1,dat.dtype);;
    dat.v_tmp = zeros(settings.N,1,dat.dtype);

    dat.i0 = trapz(dat.x,dat.J_TL);
    %current
    %     dat.i0 = dat.v0/dat.Z0;
    %current modulation function
    dat.i0_t = @(tm) dat.i0*(1+settings.modA*sin(2*pi*settings.modF*dat.f_R*tm));
    dat.i_TL = dat.i0*(settings.Ltot-dat.x)/settings.Ltot;
%     dat.i_TL = dat.i0*ones(settings.N,1);
    dat.Acoeff = (settings.nTHz/settings.nRF)^2; dat.Bcoeff = dat.dt^2/eps_ch/(settings.nRF^2);
    
    dat.Dcoeff = (settings.nTHz/settings.nRF)/dat.Z0;
    dat.Ecoeff = (settings.nTHz/settings.nRF)*dat.Z0;
    dat.Fcoeff = dat.dt/eps_ch/(settings.nRF^2);

    rho_Au = 2.3*1E-8; Ln = 300*1E-6; Lp = settings.Lp*1E-9; mu0 = Constants('mu0');
    
    dat.R_Au = rho_Au/(Ln*Lp*mu0); %in units of 1/s
    dat.R_Au = 0*dat.R_Au/1E12; 
    
end
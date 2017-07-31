function [ dat ] = makeTransLineVarsx2(params,dat)

    %convert the curr density in A/m^2 into A/mm^2 by dividing by 1E6;
    dat.J_TL1 =  zeros(params.N,1);
    dat.J_TL2 =  zeros(params.N,1);
    dat.eps_ch = 1E15*Constants('eps0'); 
    mu_ch  = 1E3*Constants('mu0');

    %   TL impedence
    dat.Z0 = sqrt(mu_ch/dat.eps_ch)/params.nRF;
    dat.width_mm = 20e-3; % width of the laser in mm
    dat.Lp_mm = params.Lp*1e-6;%
    dat.Np = 50; %nr of periods!
    dat.height_mm = dat.Lp_mm*dat.Np;


    dat.Z_kV_A = 50*1E-3; % kV/A <- connecting wire impedance.
    dat.V_s_0 = params.voltage*1e-3; % V->kV
    dat.Vs = @(tm) dat.V_s_0*(1+params.modA*sin(2*pi*tm*params.f_R*params.modF)); %12volts?

    %   bias
    dat.v0 = params.bias/10;
    dat.v_TL = dat.v0*ones(params.N,1);
    dat.v_TL_old = dat.v_TL; dat.v_tmp = zeros(params.N+1,1);

    %   current with additional ghost cell
    dat.i0 = (dat.Vs(0)-dat.v0*dat.height_mm)/(dat.Z_kV_A)/dat.width_mm;
    xI = linspace(0,params.Ltot+params.dx,params.N-1).';
    dat.i_TL = dat.i0*(xI(end)-xI)./xI(end);

    dat.Acoeff = (params.nTHz/params.nRF)^2; 
    dat.Bcoeff = params.dt^2/dat.eps_ch/(params.nRF^2);
    dat.Dcoeff = (params.nTHz/params.nRF)/dat.Z0;
    dat.Ecoeff = (params.nTHz/params.nRF)*dat.Z0;
    dat.Fcoeff = params.dt/dat.eps_ch/(params.nRF^2);

    % resistance losses!
    rho_Au = 2.3*1E-8; Ln = 300*1E-6; Lp = params.Lp*1E-9; mu0 = Constants('mu0');
    dat.R_Au = rho_Au/(Ln*Lp*mu0); %in units of 1/s
    dat.R_Au = dat.R_Au/1E12;
    dat.J_TL = dat.v_TL*0;

end


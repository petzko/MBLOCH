function dat = makeTransLineVars(settings,dat)




    %make transline vars

    %convert the curr density in A/m^2 into A/mm^2 by dividing by 1E6;
    dat.J_TL1 =  zeros(settings.N,1,dat.dtype); 
    dat.J_TL2 =  zeros(settings.N,1,dat.dtype);
    
    eps_ch = 1E15*Constants('eps0'); mu_ch  = 1E3*Constants('mu0');

    
%%%  current, bias and voltage modulation ftions
%     FWHM = 0.5; %psec
%     stdev = FWHM/(2*sqrt(2*log(2))); % std. dev
%     gauss = @(tm) exp(-(tm-1).^2/(2*stdev^2))/(stdev*sqrt(2*pi));
        
%     dat.i0_t = @(tm) dat.i0*(1+settings.modA*sin(2*pi*settings.modF*dat.f_R*tm));
%     dat.v0_t = @(tm) dat.v0*(1+settings.modA*sin(2*pi*settings.modF*dat.f_R*tm));

%width of the transmission line is 20 micrometers -> 20e-3 millimeters;
% TL impedence 
    dat.Z0 = sqrt(mu_ch/eps_ch)/settings.nRF;
    
%   voltage
    dat.width_mm = 20e-3; % width of the laser in mm
    dat.height_mm = 10e-3; % height of the laser in mm 
    dat.Lp_mm = settings.Lp*1e-6;%
    dat.Np = 50; %nr of periods! 

    dat.Z_kV_A = 50*1E-3; % kV/A <- connecting wire impedance.  
    dat.V_s_0 = settings.voltage*1e-3; % V->kV
    dat.Vs = @(tm) dat.V_s_0*(1+settings.modA*sin(2*pi*tm*dat.f_R*settings.modF)); %12volts?
    
%   bias
    dat.v0 = settings.bias/10;
    dat.v_TL = dat.v0*ones(settings.N,1,dat.dtype); 
    dat.v_TL_old = dat.v_TL; dat.v_tmp = zeros(settings.N,1,dat.dtype);


% % %     current with additional ghost cell
    dat.i0 = (dat.Vs(0)-dat.v0*dat.Lp_mm*dat.Np)/(dat.Z_kV_A+dat.Z0)/dat.width_mm; 
    xI = linspace(0,settings.Ltot+dat.dx,settings.N+1).';
    dat.i_TL = dat.i0*(xI(end)-xI)./xI(end);

%     dat.i0 = (dat.Vs(0)-dat.v0*dat.Lp_mm*dat.Np)/(dat.Z_kV_A+dat.Z0)/dat.width_mm; 
%     dat.i_TL = dat.i0*(dat.x(end)-dat.x)./dat.x(end);


    dat.Acoeff = (settings.nTHz/settings.nRF)^2; dat.Bcoeff = dat.dt^2/eps_ch/(settings.nRF^2);
    dat.Dcoeff = (settings.nTHz/settings.nRF)/dat.Z0;
    dat.Ecoeff = (settings.nTHz/settings.nRF)*dat.Z0;
    dat.Fcoeff = dat.dt/eps_ch/(settings.nRF^2);

    % resistance losses! 
    rho_Au = 2.3*1E-8; Ln = 300*1E-6; Lp = settings.Lp*1E-9; mu0 = Constants('mu0');
    dat.R_Au = rho_Au/(Ln*Lp*mu0); %in units of 1/s
    dat.R_Au = dat.R_Au/1E12; 
    
end
clearclear; close all; clc;
ratesfile = ['fitted_data.mat']; simDataFile = [ 'OPTICAMAIN.sim']; scenarioFile = [ 'IVL.set'];
init = +1;

%parse all input files and load the scatterin rates file !
settings = parseSimParams(scenarioFile);
settings = parseSimData(simDataFile,settings);
load(ratesfile);



Voltages = [10.5:2.:12]; %in kV/cm;

Is = zeros(length(Voltages),1);
Biases = zeros(length(Voltages),1);
Volts = zeros(length(Voltages),1);
Js = zeros(length(Voltages),1);
Ls = zeros(length(Voltages),1);
IVJcounter = 1;
init_bias = 8.0/10;
ratesbias = 8.0/10;
%%
for Voltage0 = Voltages
    %%
    clc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZE DEVICE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%
    %phase velocity inside the medium ( in mm per picosecond ... )RT time and
    %frequency
    dat.c = Constants('c',{'time',settings.tch},{'length',settings.lch})/settings.nTHz; dat.T_R = 2*settings.Ltot/dat.c; dat.f_R = 1/dat.T_R;
    dat.dtype = 'double';
    %%%%dipole mtx elements (in Cnm)
    % hbar in eV-ps
    dat.hbar = Constants('hbar',{'time',settings.tch})/Constants('q0');
    dat.INJ = 1; dat.ULL = 2;dat.LLL = 3; dat.RES = 4; dat.DEPOP = 5;
    resonance = (E_fit{dat.INJ}(init_bias)-E_fit{dat.ULL}(init_bias))*1E3;
    %central field frequency.
    dat.E0 =  (E_fit{dat.ULL}(init_bias)-E_fit{dat.LLL}(init_bias))/dat.hbar;
    %cavity loss l_0 in (cm^-1) --> l_0*100 in (m^-1) --> 1 mm^-1
    %make the loss large in order to supress lasing!
    dat.l_0 = 2*settings.loss*100/(1/settings.lch);
    
    %%%%%%%%%%%%%%%%% simuation parameters %%%%%%%%%%%%%%%%%
    %grid size in x direction
    settings.N = 3000;
    dat.x = linspace(0,settings.Ltot,settings.N)';
    if strcmp(dat.dtype,'single')
        dat.x = single(dat.x);
    end
    
    dat.dx = dat.x(2) - dat.x(1); dat.dt = dat.dx/dat.c;
    dat.diffusion = 4*dat.E0/dat.c^2*settings.D*10^2/(1/settings.tch);
    
    dat.zUL = zUL_fit(init_bias);
    dat.Ncarriers = settings.dN*(100^3)*settings.Ld/settings.Lp; % cm^-3 --> m^-3; carrier density
    dat.trace_rho = ((dat.E0*1E12*dat.Ncarriers*settings.Overlap*((dat.zUL*1E-9*Constants('q0'))^2))/(Constants('eps0')*settings.nTHz*Constants('c')*Constants('hbar')))/(1/(settings.lch*settings.tch));
    
    dat.NLVLS = 5;
    dat.W = zeros(settings.N,dat.NLVLS,dat.NLVLS,dat.dtype);
    dat.G = zeros(settings.N,dat.NLVLS,dat.dtype);
    dat.dipR = ones(settings.N,1,dat.dtype);
    
    
    
    
    
    
    
    
    dat = makeMaxwellVars(settings,dat);
    %%% start makeTransLineVars
    %convert the curr density in A/m^2 into A/mm^2 by dividing by 1E6;
    dat.J_TL1 =  zeros(settings.N,1,dat.dtype);
    dat.J_TL2 =  zeros(settings.N,1,dat.dtype);
    
    eps_ch = 1E15*Constants('eps0'); mu_ch  = 1E3*Constants('mu0');
    dat.Z0 = sqrt(mu_ch/eps_ch)/settings.nRF;
    
    %   voltage
    dat.width_mm = 20e-3; % width of the laser!
    dat.Lp_mm = settings.Lp*1e-6;%
    dat.Z_kV_A = 50*1E-3; % kV/A <- connecting wire.
    dat.V_s_0 = Voltage0*1e-3; %12 V->kV;
    dat.Vs = @(tm) dat.V_s_0*(1+settings.modA*sin(2*pi*tm*dat.f_R*settings.modF)); %12volts?
    
    
    %   bias
    dat.v0 = init_bias;
    dat.v_TL = dat.v0*ones(settings.N,1,dat.dtype);
    dat.v_TL_old = dat.v_TL; dat.v_tmp = zeros(settings.N,1,dat.dtype);
    
    
    %     current
    dat.i0 = (dat.Vs(0)-dat.v0*dat.Lp_mm)/dat.Z_kV_A/dat.width_mm;
    xI = linspace(0,settings.Ltot+dat.dx,settings.N+1).';
    dat.i_TL = dat.i0*(xI(end)-xI)./xI(end);
    
    dat.Acoeff = (settings.nTHz/settings.nRF)^2; dat.Bcoeff = dat.dt^2/eps_ch/(settings.nRF^2);
    dat.Dcoeff = (settings.nTHz/settings.nRF)/dat.Z0;
    dat.Ecoeff = (settings.nTHz/settings.nRF)*dat.Z0;
    dat.Fcoeff = dat.dt/eps_ch/(settings.nRF^2);
    
    % resistance losses!
    rho_Au = 2.3*1E-8; Ln = 300*1E-6; Lp = settings.Lp*1E-9; mu0 = Constants('mu0');
    dat.R_Au = rho_Au/(Ln*Lp*mu0); %in units of 1/s
    dat.R_Au = dat.R_Au/1E12;
    %%% end makeTransLineVars
    
    
    
    dat.t = dat.dt;
    dat.iter_ctr = 0;
    
    dat = interpParams(settings,dat,W_fit,E_fit,AC_fit,zUL_fit);
    dat = makeBlochVars(settings,dat);
    
    dat.t = dat.dt;
    
    
    INJ = dat.INJ; ULL = dat.ULL; LLL = dat.LLL; RES = dat.RES; DEPOP = dat.DEPOP;
  
    
    
    
    
    %this is needed for current density calculations!
    %     dat.rates =   dat.r110.*dat.W(:,dat.INJ,dat.DEPOP) + dat.r220.*dat.W(:,dat.LLL,dat.DEPOP) + ...
    %     dat.r330.*dat.W(:,dat.ULL,dat.DEPOP)+dat.rRES.*dat.W(:,dat.RES,dat.DEPOP);
    %     dat.rates_t = zeros(settings.N,1,dat.dtype);
    
    
    
    
    %%%% specify some of the mainloop control parameters %%%%
    iter_per_rt = round(dat.T_R/dat.dt);
    
    plotCtr = 500;
    interpCtr = 100;
    maxITER = 5E6;
    means = zeros(maxITER,1);
    mean_ctr = 1;
    settings.simRT = maxITER*dat.dt/dat.T_R;
    %simulation info storage arrays -> preallocate
    
    info.settings = settings;
    info.cavity = 'FP-OPTICA';
    info.Ltot = settings.Ltot;
    info.N = settings.N;
    info.SIMTYPE = ['TL with (i0,v0) = (' num2str(dat.i0*10) ',' num2str(init_bias*10) ') in units (A/cm,kV/cm)'];
    unfinished = true;
    v_TL_old = zeros(settings.N,1);
    
    while (dat.iter_ctr < maxITER && unfinished)
        
        if(mod(dat.iter_ctr+1,interpCtr) == 0)
            dat = interpParams(settings,dat,W_fit,E_fit,AC_fit,zUL_fit);
        end
        
        if (mod(dat.iter_ctr+1,100) == 0)
            
            m_n = mean(dat.v_TL) ;
            means(mean_ctr) = m_n;
            mean_ctr = mean_ctr+1;
            if abs(mean(v_TL_old)-m_n)< 1E-8
                Biases(IVJcounter) = m_n;
                Volts(IVJcounter) = Voltage0;
                Is(IVJcounter) = mean(dat.i_TL);
                Js(IVJcounter) = mean(dat.J_TL);
                Ls(IVJcounter) = mean(abs(dat.U).^2+abs(dat.V).^2);
                IVJcounter = IVJcounter+1;
                unfinished = false;
            end
            
            if(m_n < 6/10 || m_n>13/10)
                unfinished = false;
                IVJcounter = IVJcounter+1;
            end
            v_TL_old =dat.v_TL;
            %set the initial voltage to be the current mean value of the
            %voltage. should speed up convergence!
            init_bias = m_n;
        end
        
        %%plot some of the results if neeed ariseth :D
        if(mod(dat.iter_ctr,100) == 0)
            clc;
            info.iter_ctr = dat.iter_ctr;
            info.RT = dat.t/dat.T_R;
            intensity = dat.U.*conj(dat.U) + dat.V.*conj(dat.V) ;
            info.maxInt  =  max(intensity);
            printINFO(info);
            clc;
            info.dat.iter_ctr = dat.iter_ctr;
            info.RT = dat.t/dat.T_R;
            intensity = dat.U.*conj(dat.U) + dat.V.*conj(dat.V) ;
            info.maxInt  =  max(intensity);
            printINFO(info);
            
            display(['init dE_{1''3}: ' num2str(resonance)]);
            display(['V_s: ' num2str(dat.Vs(dat.t))]);
            display(['v_0: ' num2str(dat.v_TL(1))]);
            display(['i-in: ' num2str(dat.i_TL(1))]);
            display(['i-out1: ' num2str(trapz(dat.x,dat.J_TL1))]);
            display(['i-out2: ' num2str(trapz(dat.x,dat.J_TL2))]);
            
            subplot(3,1,1)
            plot(dat.x,[real(dat.U),real(dat.V)]);
            subplot(3,1,2)
            plotyy(dat.x,[dat.v_TL*10],[dat.x;dat.x(end)+dat.dx],dat.i_TL);
            title(info.SIMTYPE);
            subplot(3,1,3)
            %plots the populations and the current density
            plot(dat.x,dat.J_TL1,dat.x,dat.J_TL2);
            ylim([1,5]);
            getframe;
        end
        
        dat = stepBloch(settings,dat);
        dat = stepTransLine(settings,dat);
        dat = stepWave(settings,dat);
        dat = updateBloch(settings,dat);
        
        dat.t = dat.t+dat.dt;
        dat.iter_ctr = dat.iter_ctr + 1;
    end
    means= means(1:mean_ctr-1);
    save(['means_95REFL' num2str(Voltage0)],'means','Voltage0');
    
end
save('IVCHARS_FINAL');
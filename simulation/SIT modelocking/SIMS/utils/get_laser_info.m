function [ dat ] = get_laser_info( sim_settings )

    %length of the absorption and the gain region
    dat.Ltot =  sim_settings.Lg + sim_settings.La;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZE DEVICE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%
    % central freq. phase velocity in characteristic units
    dat.c = Constants('c',{'time',sim_settings.tch},{'length',sim_settings.lch})/sim_settings.nTHz;
    % round trip time and frequency
    dat.T_R = dat.Ltot/dat.c; dat.f_R = 1/dat.T_R;
    % hbar in eV-ps
    hbar = Constants('hbar',{'time',sim_settings.tch})/Constants('q0');
    
    
    %indices of the upper and lower laser level
    ULL = sim_settings.ULL ; LLL = sim_settings.LLL;
    %indices of the gain and absorptopn secions
    GAIN = 1; ABS = 2;
    %eigenenergies of the gain hamiltonian -> transform into ang. freq.
    NLVLS = 2;
    Hg = reshape(sim_settings.Hg',NLVLS,NLVLS);
    E2g = Hg(ULL,ULL)/hbar;
    E1g = Hg(LLL,LLL)/hbar;
    
    %eigenenergies of the gain hamiltonian -> transform into ang. freq.
    Ha = reshape(sim_settings.Ha',NLVLS,NLVLS);
    E2a = Ha(ULL,ULL)/hbar;
    E1a = Ha(LLL,LLL)/hbar;
    
    E21g = E2g-E1g; %rad/ps  2->1 gain transition freq (optical transition)
    E21a = E2a-E1a; %rad/ps; 2->1 absorber transition freq
    %central field frequency.
    dat.E0 = (E21g+E21a)/2;
    
    %cavity loss l_0 in (cm^-1) --> l_0*100 in (m^-1) --> 1 mm^-1
    dat.lg = sim_settings.gain_loss*100/(1/sim_settings.lch);
    dat.la = sim_settings.abs_loss*100/(1/sim_settings.lch);
    
    %gain section carrier density (doping density) and trace normalization
    %constant
    dat.Ncarriers_g = sim_settings.dN_g*(100^3)*sim_settings.Ld_g/sim_settings.Lp_g; % cm^-3 --> m^-3; carrier density
    %absorption section carrier density and trace normalization constant
    dat.Ncarriers_a = sim_settings.dN_a*(100^3)*sim_settings.Ld_a/sim_settings.Lp_a; % cm^-3 --> m^-3; carrier density
       
   mu_g = sim_settings.zULg*1E-9*Constants('q0');
    mu_a = sim_settings.zULa*1E-9*Constants('q0');

    sigma_g = sim_settings.T2_g*(1e-12)*sim_settings.Overlap_g*dat.E0*1E12*mu_g^2/ ...
        (Constants('eps0')*sim_settings.nTHz*Constants('c')*Constants('hbar'));
    sigma_a = sim_settings.T2_a*(1e-12)*sim_settings.Overlap_a*dat.E0*1E12*mu_a^2/ ...
        (Constants('eps0')*sim_settings.nTHz*Constants('c')*Constants('hbar'));

    Ga = sigma_a*dat.Ncarriers_a*sim_settings.La*1e-3;
    Gg = sigma_g*dat.Ncarriers_g*sim_settings.Lg*1e-3;
    powerloss = 2*sim_settings.gain_loss*100;
    dat.d_th =  powerloss/(sigma_g*dat.Ncarriers_g)+Ga/Gg;
    
    dat.NLVLS = 2;
    
    %%%
    
    
    invth = dat.d_th; 
    dat.thgain = (dat.E0*sim_settings.T2_g/(Constants('hbar')*Constants('eps0')*...
        Constants('c')*sim_settings.nTHz))*dat.Ncarriers_g*...
        (sim_settings.zULg*1e-9*Constants('q0'))^2*invth*1e-2; % per cm
    
    invmax = 1;
    dat.maxgain = (dat.E0*sim_settings.T2_g/(Constants('hbar')*Constants('eps0')*...
        Constants('c')*sim_settings.nTHz))*dat.Ncarriers_g*...
        (sim_settings.zULg*1e-9*Constants('q0'))^2*invmax*1e-2; % per cm

end


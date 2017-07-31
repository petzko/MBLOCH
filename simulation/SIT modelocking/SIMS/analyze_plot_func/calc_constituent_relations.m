
% this is in order to calculate the output power
Width = 40*1E-6; % 20um 
Height = 10e-6; % 10um 

% The e-field is normalized 
zUL = sim_settings.zULg; 
nthz = sim_settings.nTHz;
normalization_factor = 1E12*Constants('hbar')/(Constants('q0')*zUL*1E-9);

dt = dat.dt; 
iter_per_rt = round(dat.T_R/dat.dt);

envelope = record_U_a;

lbdry = 290*iter_per_rt; rbdry = 298*iter_per_rt;
intrv1 = [lbdry:rbdry];

p = sim_settings.r22g_0/d_th;

mu_g = sim_settings.zULg*1E-9*Constants('q0');
mu_a = sim_settings.zULa*1E-9*Constants('q0');

sigma_g = sim_settings.T2_g*(1e-12)*sim_settings.Overlap_g*dat.E0*1E12*mu_g^2/ ...
    (Constants('eps0')*sim_settings.nTHz*Constants('c')*Constants('hbar'));
sigma_a = sim_settings.T2_a*(1e-12)*sim_settings.Overlap_a*dat.E0*1E12*mu_a^2/ ...
    (Constants('eps0')*sim_settings.nTHz*Constants('c')*Constants('hbar'));


% this is the simulation data intensity normalized to 1 for max intensity
data1 = abs(envelope(intrv1)).^2;
data1 = data1/max(data1);
% find the peaks and the FWHM values! 
[pks,locs,w,heights] = findpeaks(data1,'MinPeakHeight',0.3,'MinPeakProminence',0.5);


I0 = max(abs(envelope(intrv1))).^2;
I0_ = Constants('eps0')*sim_settings.nTHz*Constants('c')/8*(abs((I0)*normalization_factor).^2);
results.I0 = I0;
results.T1g = sim_settings.T1_g;
results.T2g = sim_settings.T2_g;
results.p = p; 
results.d_th = d_th;
results.P0_mW = I0_*Width*Height*1e3; results.p = sim_settings.r22g_0/d_th; results.d_th = d_th; 


% remaining constituent eqns 
results.gth = d_th*sigma_a*dat.Ncarriers_a*sim_settings.La*1e-3; % unitless
results.g0 = p*results.gth; % unitless
results.q0 =  results.p*results.d_th*sigma_g*dat.Ncarriers_g*sim_settings.La*1e-3; %unitless

results.dipole_ratio = sim_settings.zUL_g/sim_settings.zUL_a;

WGsat = 1./(sim_settings.T1_g*sim_settings.T2_g);
WAsat = results.dip_R^2./(sim_settings.T1_a*sim_settings.T2_a);

display(['q0/WAsat-g0/WGsat:', num2str(results.q0/(sim_settings.T1_a*WAsat)^2 - results.g0/(sim_settings.T1_g*WGsat)^2)]);

kappa_a = 4;
kappa_g = 2; 

display(['q0*T2a^2*kappaa^4-g0*T2g^2*kappag^4: ' , num2str(results.q0*sim_settings.T2_a^2*kappa_a^4-results.g0*sim_settings.T2_g^2*kappa_g^4)]);
results.aLg = 2*sim_settings.gain_loss*100*sim_settings.Lg*1e-3; % unitless
results.kappaA = mu_a/Constants('hbar'); %m/Vs%
results.kappaG = mu_g/Constants('hbar'); %m/Vs%


results.WA_sat = 1/((results.kappaA)^2*sim_settings.T1_a*sim_settings.T2_a*1e-24); %V^2/m^2
results.WG_sat = 1/((results.kappaG)^2*sim_settings.T1_g*sim_settings.T2_g*1e-24); %V^2/m^2
results.gamma = results.q0/results.WA_sat; 


if length(pks) == 0
    results.pulsed = false;
    display('analysis done no pulsations');
    return 
end

results.higher_harmonic = (locs(2)-locs(1))*dt < 0.75*dat.T_R;

results.sim_settings = sim_settings;
results.dt = dat.dt;
results.pulsed = true;

results.FWHM = w(2)*dat.dt;  % FWHM of the simulated pulse 
results.tau_real = results.FWHM/1.763; % tau of the simulated pulse see A. Weiners "Ultrafast optics"
results.deltaT = sim_settings.T1_g*log((p*d_th -...
    (record_r22g(intrv1) - record_r11g(intrv1)))./(d_th*(p-1)));

WA_sat = 1/(sim_settings.T1_a*sim_settings.T2_a); % 1/ps^2
gamma = results.q0/WA_sat; % ps^2

results.tau_est = sim_settings.T2_g./sqrt((results.p-1)*results.gth); % ps
I0_est1 = sqrt((results.p-1)*results.gth)/sim_settings.T2_g.^2; %1/ps^2
I0_est2 = 1/4*(5*results.p-4)*results.gth/gamma;  %1/ps^2

I0_est_1 = Constants('eps0')*sim_settings.nTHz*Constants('c')/8*(abs((I0_est1)*normalization_factor).^2);
I0_est_2 = Constants('eps0')*sim_settings.nTHz*Constants('c')/8*(abs((I0_est2)*normalization_factor).^2);
   
% options = optimoptions('fsolve','Display','final-detailed','PlotFcn',@optimplotfirstorderopt);
% options.MaxFunEvals = 10000;
% options.MaxIter = 10000;
options = [];

params.T1_g = sim_settings.T1_g; %ps 
params.T2_g = sim_settings.T2_g; %ps
% dipole ratio 
r = sim_settings.zULg/sim_settings.zULa;
params.WG_sat = 1/(sim_settings.T1_g*sim_settings.T2_g); %ps^-2;
params.WA_sat = r^2/(sim_settings.T1_a*sim_settings.T2_a); %ps^-2;

params.g0 = results.g0;
params.gamma = gamma;
params.q0 = results.q0;

params.p = p; 
params.d_th = d_th;

[estimates,fval] = fsolve('constituent_relations_FSA_SG',[0,I0,results.tau_real],options,params); 
display(fval)
results.estimates = estimates;
I0_est_3 = Constants('eps0')*sim_settings.nTHz*Constants('c')/8*(abs((estimates(2))*normalization_factor).^2);

results.P0_est_mW1 = I0_est_1*Width*Height*1e3; %mW
results.P0_est_mW2 = I0_est_2*Width*Height*1e3; %mW
results.P0_est_mW3 = I0_est_3*Width*Height*1e3; %mW

display('analysis done');

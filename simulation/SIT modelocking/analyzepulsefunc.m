function [ results ] = analyzepulsefunc( datafile )
    
load(datafile);

% typical QCL gain medium widht and height in transverse direction -> 
% width x height is the outcoupling area used to calculate the outcoupled 
% power! 

Width = 40*1E-6; % 20um 
Height = 10e-6; 
% the electric field is normalized to the dipole / hbar moment of the gain
% material. Thus now the field has units of freq. (the rabi freq). 
%!!! save this normalization constant 
zUL = sim_settings.zULg;

% this number is used to undo the normalization 
normalization_factor = 1E12*Constants('hbar')/(Constants('q0')*zUL*1E-9);

% time step (ps)
dt = dat.dt; 
% iterations per round trip (roughly) 
iter_per_rt = round(dat.T_R/dat.dt);

% this is the vector holding the rabi-field envelope
envelope = record_U_a;

% sub-set of the envelope -> take the last couple of roudn trips only since
% only then will the pulse have converged 
lbdry = 290*iter_per_rt; rbdry = 298*iter_per_rt;
intrv1 = [lbdry:rbdry];

% pump power. d_th is the threshold inversion (calculated via Eq. (17)). 
p = sim_settings.r22g_0/d_th;

% tipole moments in SI units for the gain (g) and absorber (a) 
mu_g = sim_settings.zULg*1E-9*Constants('q0');
mu_a = sim_settings.zULa*1E-9*Constants('q0');

% gain (g) and absorber cross-section for the central frequency of the
% pulse. 
sigma_g = sim_settings.T2_g*(1e-12)*sim_settings.Overlap_g*dat.E0*1E12*mu_g^2/ ...
    (Constants('eps0')*sim_settings.nTHz*Constants('c')*Constants('hbar'));
sigma_a = sim_settings.T2_a*(1e-12)*sim_settings.Overlap_a*dat.E0*1E12*mu_a^2/ ...
    (Constants('eps0')*sim_settings.nTHz*Constants('c')*Constants('hbar'));


% calculate the intensity (in SI units) of the pulse. 
I0 = max(abs(envelope(intrv1))).^2;
I0_ = Constants('eps0')*sim_settings.nTHz*Constants('c')/8*(abs((I0)*normalization_factor).^2);

results.I0 = I0;
results.T1g = sim_settings.T1_g;
results.T2g = sim_settings.T2_g;
results.p = p; 
results.d_th = d_th;
results.P0_mW = I0_*Width*Height*1e3; results.p = sim_settings.r22g_0/d_th; results.d_th = d_th; 

% threshold gain x propagation length
results.gth = d_th*sigma_g*dat.Ncarriers_g*sim_settings.Lg*1e-3; % unitless
% small signal gain  
results.g0 = p*results.gth; % unitless
% small signal losses (absorber)
results.q0 =  sigma_a*dat.Ncarriers_a*sim_settings.La*1e-3; %unitless

% ratio of the dipole moments 
results.dipole_ratio = sim_settings.zULg/sim_settings.zULa;
% saturation intensity (normalized to rabi freq. ^2) of the gain and the 
% absorber
WGsat = 1./(sim_settings.T1_g*sim_settings.T2_g);
WAsat = results.dipole_ratio^2./(sim_settings.T1_a*sim_settings.T2_a);

display(['q0/T1a*WAsat-g0/T1g*WGsat:', num2str(results.q0/(sim_settings.T1_a*WAsat)^2 - results.g0/(sim_settings.T1_g*WGsat)^2)]);


kappa_a = 4;
kappa_g = 2; 

display(['q0*T2a^2*kappaa^4-g0*T2g^2*kappag^4: ' , num2str(results.q0*sim_settings.T2_a^2*kappa_a^4-results.g0*sim_settings.T2_g^2*kappa_g^4)]);
results.aLg = 2*sim_settings.gain_loss*100*sim_settings.Lg*1e-3; % unitless
results.kappaA = mu_a/Constants('hbar'); %m/Vs% % coupling factor absorber
results.kappaG = mu_g/Constants('hbar'); %m/Vs% % coupling factor gain 

% saturation e-field^2 in SI units 
results.WA_sat = 1/((results.kappaA)^2*sim_settings.T1_a*sim_settings.T2_a*1e-24); %V^2/m^2
results.WG_sat = 1/((results.kappaG)^2*sim_settings.T1_g*sim_settings.T2_g*1e-24); %V^2/m^2
% self amplitude modulation coefficient. 
results.gamma = results.q0/results.WA_sat; 

%%%% NOW ANALYZE THE SIMULATION DATA (E-field envelope)

% this is the simulated envelope's intensity normalized to 1 for max intensity
data1 = abs(envelope(intrv1)).^2;
data1 = data1/max(data1);

% find the peaks and the FWHM values! 
[pks,locs,w,heights] = findpeaks(data1,'MinPeakHeight',0.3,'MinPeakProminence',0.5);

% check if there were pulses!
if length(pks) == 0
    results.pulsed = false;
    display('analysis done no pulsations');
    return 
end

% check if the separation of two consecutive pulses is less than the 75% of
% the round trip time (second and higher harmonic mode-locking) 
results.higher_harmonic = (locs(2)-locs(1))*dt < 0.75*dat.T_R;

results.sim_settings = sim_settings;
results.dt = dat.dt;
results.pulsed = true;

% FWHM of the simulated pulse 
results.FWHM = w(2)*dat.dt;  
% tau of the simulated pulse assuming a sech pulse! that  see A. Weiners "Ultrafast optics"
results.tau_real = results.FWHM/1.763; 
% gain recovery time ( see Eq. (21) ) ! 
results.deltaT = sim_settings.T1_g*log((p*d_th -...
    (record_r22g(intrv1) - record_r11g(intrv1)))./(d_th*(p-1)));

WA_sat = results.dipole_ratio/(sim_settings.T1_a*sim_settings.T2_a); % 1/ps^2
gamma = results.q0/WA_sat; % ps^2 -> recalculate gamma in normalized units! 

% ps -> estimate of tau in case alpha = 0 ( from Eq. (36) ) 
results.tau_est = sim_settings.T2_g./sqrt((results.p-1)*results.gth);
% estimate of I0 in case alpha = 0
I0_est_1 = sqrt((results.p-1)*results.gth)/sim_settings.T2_g.^2; %1/ps^2
I0_est_2 = 1/4*(5*results.p-4)*results.gth/gamma;  % in units of 1/ps^2

I0_est_1 = Constants('eps0')*sim_settings.nTHz*Constants('c')/8*(abs((I0_est_1)*normalization_factor).^2);
I0_est_2 = Constants('eps0')*sim_settings.nTHz*Constants('c')/8*(abs((I0_est_2)*normalization_factor).^2);

%%%%% Calculate the estimates for E0^2 and tau from Eqs. (36)-(38).

% use Matlab's fsolve to treat the more general case when alpha neq 0. 

% % custom options for the solver (in case convergence is not achieved) 
% options = optimoptions('fsolve','Display','final-detailed','PlotFcn',@optimplotfirstorderopt);
options = optimoptions('fsolve');
options.MaxFunEvals = 10000;
options.MaxIter = 10000;

% default options for the solver
% options = [];

% some params for the constituent relations calculator
params.T1_g = sim_settings.T1_g; %ps 
params.T2_g = sim_settings.T2_g; %ps
params.WG_sat = 1/(sim_settings.T1_g*sim_settings.T2_g); %ps^-2;
params.WA_sat = results.dipole_ratio^2/(sim_settings.T1_a*sim_settings.T2_a); %ps^-2;

params.g0 = results.g0;
params.gamma = gamma;
params.q0 = results.q0;

params.p = p; 
params.d_th = d_th;

[estimates,fval] = fsolve('constituent_relations_FSA_SG',[0,0,results.tau_real],options,params); 
display(fval)

results.estimates = estimates;

I0_est_3 = Constants('eps0')*sim_settings.nTHz*Constants('c')/8*(abs((estimates(2))*normalization_factor).^2);
results.P0_est_mW1 = I0_est_1*Width*Height*1e3; %mW
results.P0_est_mW2 = I0_est_2*Width*Height*1e3; %mW
results.P0_est_mW3 = I0_est_3*Width*Height*1e3; %mW

display('analysis done');

end


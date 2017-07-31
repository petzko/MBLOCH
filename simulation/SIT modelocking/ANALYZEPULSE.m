close all; clear;
%%%% PULSE ANALYSIS SCRIPT 
% load data
% - extract simulated pulse 
% - extract FWHM 
% - calculate gain recovery time 
% - calculate the interferometric autocorrelation trace of the pulse and
% compare it to the IAC trace of a gaussian or a sech pulse with the same
% FWHM! 

load('SIMRES/test1-mod1p1_N_4000_300.mat')


dt = dat.dt;
envelope = record_U_a;


iter_per_rt = round(dat.T_R/dat.dt);

lbdry = 290*iter_per_rt;
rbdry = 298*iter_per_rt;
intrv1 = [lbdry:rbdry];
% this is the simulation data intensity normalized to 1 for max intensity
data1 = abs(envelope(intrv1)).^2;
data1 = data1/max(data1);
% times vector
tmsdata1 = [-length(data1)/2:length(data1)/2-1].'*dat.dt;
% plot(tmsdata1/dat.T_R,data1);
% find the peaks and the FWHM values! 
[pks,locs,w,heights] = findpeaks(data1,'MinPeakHeight',0.01);
FWHM = w(3)*dat.dt;  % FWHM of the simulated pulse 
tau = FWHM/1.763; % tau of the simulated pulse see A. Weiners "Ultrafast optics"
loc_I = locs(4); 
dLoc = locs(4)-locs(3); 

higher_harmonic = (locs(2)-locs(1))*dt < 0.75*dat.T_R;



% get only the middle pulse 
interv2 = [loc_I-round(dLoc/2):loc_I+round(dLoc/2)]
data2 = data1(interv2);
data2 = data2/max(data2); 

tmsdata2 = [-length(data2)/2:length(data2)/2-1].'*dat.dt;


I0 = max(abs(envelope(intrv1))).^2;

Wsat = 1/(sim_settings.T1_g*sim_settings.T2_g);
LHS = I0;  RHS = sim_settings.T1_g*Wsat/tau;
display(['I0 = ' num2str(LHS), '; T1/tau*Wsat = ' num2str(RHS) ])
display(['RHS/LHS = ', num2str(RHS/LHS)]); 


testtms = linspace(tmsdata2(1),tmsdata2(end),length(data2));

FWHM_sech = FWHM;
tp = FWHM_sech/1.763;
test_sech_pulse = sech(testtms/tp).^2;

% make a gaussian pulse with the same FWHM
FWMH_gauss = FWHM;
tp = FWMH_gauss/1.177;
test_gauss_pulse = exp(-testtms.^2/(tp^2)).^2;


[ AC1,times_ac1,dummy1,dummy2] = Interferometric_AC(test_sech_pulse.* ...
    exp(-1i*dat.E0*testtms),testtms(2)-testtms(1) );
[ AC2,times_ac2,dummy1,dummy2] = Interferometric_AC(test_gauss_pulse.* ...
    exp(-1i*dat.E0*testtms),testtms(2)-testtms(1) );
[ AC,times_ac,dummy1,dummy2] = Interferometric_AC(data2.* ...
    exp(-1i*dat.E0*tmsdata2),testtms(2)-testtms(1) );


figsize = [0,0,0.2,0.2]
figure('units','normalized','position',figsize)
subplot(2,1,1); 
plot(tmsdata2,data2,'-k',testtms,test_sech_pulse,'-r',testtms,test_gauss_pulse,'-b');
legend('sim. data','sech','gauss');
xlabel('delay (ps)');
ylabel('Intensity');

subplot(2,1,2); 
plot(times_ac,AC,'-k',times_ac1,AC1,'-r',times_ac2,AC2,'-b');
legend('sim. data','sech','gauss');
xlabel('delay (ps)');
ylabel('IAC trace');


delta_0 = record_r22g(intrv1) - record_r11g(intrv1) ;
% dat2 = get_laser_info( sim_settings )
% p = 1/dat2.d_th;
% d_th = dat2.d_th;

deltaT = sim_settings.T1_g*log((p*d_th-delta_0)./(d_th*(p-1)));

figsize = [0,0,0.2,0.2]
figure('units','normalized','position',figsize)
ax = plotyy(tmsdata2,deltaT(interv2),tmsdata2,data2);
set(ax(1).YLabel,'String','\tau_{gr} (ps)');
set(ax(2).YLabel,'String','Pulse intensity (a.u.)');
set(ax(1).XLabel,'String','Time (ps)');

ax(1).Children.Color = [1,0,0]
ax(2).Children.Color = [0,0,1]
ax(1).YColor = [0,0,0];
ax(2).YColor = [0,0,0];




% remaining constituent eqns 
mu_g = sim_settings.zULg*1E-9*Constants('q0');
mu_a = sim_settings.zULa*1E-9*Constants('q0');

sigma_g = sim_settings.T2_g*(1e-12)*sim_settings.Overlap_g*dat.E0*1E12*mu_g^2/ ...
    (Constants('eps0')*sim_settings.nTHz*Constants('c')*Constants('hbar'));
sigma_a = sim_settings.T2_a*(1e-12)*sim_settings.Overlap_a*dat.E0*1E12*mu_a^2/ ...
    (Constants('eps0')*sim_settings.nTHz*Constants('c')*Constants('hbar'));

Ga = sigma_a*dat.Ncarriers_a*sim_settings.La*1e-3;
Gg = sigma_g*dat.Ncarriers_g*sim_settings.Lg*1e-3;
powerloss = 2*sim_settings.gain_loss*100;
d_th =  powerloss/(sigma_g*dat.Ncarriers_g)+Ga/Gg
p = sim_settings.r22g_0/d_th

gth = d_th*sigma_a*dat.Ncarriers_a*sim_settings.La*1e-3
g0 = p*gth


tau_est = sim_settings.T2_g./sqrt((p-1)*gth);
tau_real = tau;


I0 = max(abs(envelope(intrv1))).^2;
I0_est = (1/sim_settings.T2_g^2)*sqrt((p-1)*gth);

display(['estimated I0 ' num2str(I0_est)]);
display(['actual I0 ' num2str(I0)]);
display(['act/est ratio ' num2str(I0/I0_est)])



display(['estimated tau ' num2str(tau_est)]);
display(['actual tau ' num2str(tau_real)]);
display(['act/est ratio ' num2str(tau_real/tau_est)])

q0 =  p*d_th*sigma_g*dat.Ncarriers_g*sim_settings.Lg*1e-3;
aLg = 2*sim_settings.gain_loss*100*sim_settings.Lg*1e-3;

kappaA = mu_a/Constants('hbar');
WA_sat = 1/((kappaA)^2*sim_settings.T1_a*sim_settings.T2_a*1e-24);




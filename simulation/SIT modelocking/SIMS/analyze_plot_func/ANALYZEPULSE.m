close all

% dt = dat.dt;
envelope = record_U_a;


iter_per_rt = round(T_RT/dt);

lbdry = 300*iter_per_rt;
rbdry = 308*iter_per_rt;
intrv1 = [lbdry:rbdry];
% this is the simulation data intensity normalized to 1 for max intensity
data1 = abs(envelope(intrv1)).^2;
data1 = data1/max(data1);
% times vector
tmsdata1 = [-length(data1)/2:length(data1)/2-1].'*dt;
% plot(tmsdata1/dat.T_R,data1);
% find the peaks and the FWHM values!
[pks,locs,w,heights] = findpeaks(data1,'MinPeakHeight',0.01);
FWHM = w(3)*dt;  % FWHM of the simulated pulse
tau = FWHM/1.763; % tau of the simulated pulse see A. Weiners "Ultrafast optics"
loc_I = locs(4);
dLoc = locs(4)-locs(3);

higher_harmonic = (locs(2)-locs(1))*dt < 0.75*T_RT;



% get only the middle pulse
interv2 = [loc_I-round(dLoc/2):loc_I+round(dLoc/2)]
data2 = data1(interv2);
data2 = data2/max(data2);

tmsdata2 = [-length(data2)/2:length(data2)/2-1].'*dt;


I0 = max(abs(envelope(intrv1))).^2;

Wsat = 1/(params_gain.T_1*params_gain.T_2);
LHS = I0;  RHS = params_gain.T_1*Wsat/tau;
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

E0 = (gain_model.E0+abs_model.E0)/2;
[ AC1,times_ac1,dummy1,dummy2] = Interferometric_AC(test_sech_pulse.* ...
    exp(-1i*E0*testtms),testtms(2)-testtms(1) );
[ AC2,times_ac2,dummy1,dummy2] = Interferometric_AC(test_gauss_pulse.* ...
    exp(-1i*E0*testtms),testtms(2)-testtms(1) );
[ AC,times_ac,dummy1,dummy2] = Interferometric_AC(data2.* ...
    exp(-1i*E0*tmsdata2),testtms(2)-testtms(1) );


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

deltaT = params_gain.T_1*log((pump_strength*d_th-delta_0)./(d_th*(pump_strength-1)));

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
mu_g = params_gain.zUL*1E-9*Constants('q0');
mu_a = params_abs.zUL*1E-9*Constants('q0');

sigma_g = params_gain.T_2*(1e-12)*params_gain.Gamma*E0*1E12*mu_g^2/ ...
    (Constants('eps0')*params_gain.nTHz*Constants('c')*Constants('hbar'));
sigma_a = params_abs.T_2*(1e-12)*params_abs.Gamma*E0*1E12*mu_a^2/ ...
    (Constants('eps0')*params_abs.nTHz*Constants('c')*Constants('hbar'));
Ncarriers_a = params_abs.Ncarriers_cm*1e6;
Ncarriers_g = params_gain.Ncarriers_cm*1e6;

Ga = sigma_a*Ncarriers_a*params_abs.L*1e-3;
Gg = sigma_g*Ncarriers_g*params_gain.L*1e-3;
powerloss = 2*params_gain.linear_loss*100;
d_th =  powerloss/(sigma_g*Ncarriers_g)+Ga/Gg
p = (gain_model.rho_u_0-gain_model.rho_l_0)/d_th


g0 =  p*d_th*sigma_g*Ncarriers_g*params_gain.L*1e-3;
aLg = 2*params_gain.linear_loss*100*params_gain.L*1e-3;

tau_est = params_gain.T_2./sqrt((p-1)*gth);
tau_real = tau;


I0 = max(abs(envelope(intrv1))).^2;
I0_est = (1/params_gain.T_2^2)*sqrt((p-1)*gth);

display(['estimated I0 ' num2str(I0_est)]);
display(['actual I0 ' num2str(I0)]);
display(['act/est ratio ' num2str(I0/I0_est)])



display(['estimated tau ' num2str(tau_est)]);
display(['actual tau ' num2str(tau_real)]);
display(['act/est ratio ' num2str(tau_real/tau_est)])



kappaA = mu_a/Constants('hbar');
WA_sat = 1/((kappaA)^2*params_abs.T_1*params_abs.T_2*1e-24);

% finally plot the pulse MAX/FWHM over the cavity length
iter_per_rt
I0_of_x = [];

FWHM_of_x = [];
Width = 40*1E-6; % 20um
Height = 10e-6;
if exist('pulse_data_vs_x')
    PULSE_ENERGY_of_x = zeros(size(pulse_data_vs_x,1),1);
    E0_of_x = zeros(size(pulse_data_vs_x,1),1); 
    
    for x_idx = 1:size(pulse_data_vs_x,1)
        x_idx
        envelope_x = pulse_data_vs_x(x_idx,:);
        tmsdata_x = [-length(envelope_x)/2:length(envelope_x)/2-1].'*dt;
        norm_envelope = envelope_x/max(abs(envelope));
        
        [pulse_peaks,locs] = findpeaks(abs(norm_envelope),tmsdata_x,'MinPeakHeight',0.5);
        dlocs = locs(ceil(length(locs)/2)) - locs(ceil(length(locs)/2)-1);
        peak_idx = locs(ceil(length(locs)/2))-dlocs/2 < tmsdata_x  & tmsdata_x < locs(ceil(length(locs)/2)) + dlocs/2;
        peak_env  = envelope_x(peak_idx);
        peak_tm = tmsdata_x(peak_idx);
        
        PULSE_ENERGY_of_x(x_idx) = 1./(peak_tm(end)-peak_tm(1))*trapz(peak_tm,abs(peak_env).^2);
              
        % plot(tmsdata1/dat.T_R,data1);
        % find the peaks and the FWHM values!
        
        [pks,locs,w,heights] = findpeaks(abs(envelope_x).^2/max(abs(envelope_x).^2),'MinPeakHeight',0.01);
        
        FWHM_x = max(w)*dt;  % FWHM of the simulated pulse
        tau_x = FWHM_x/1.763; % tau of the simulated pulse see A. Weiners "Ultrafast optics"
        normalization_factor = 1E12*Constants('hbar')/(Constants('q0')*params_gain.zNORM*1E-9);
        I0_x = max(abs(envelope_x(1:iter_per_rt))).^2;
        I0_x = Constants('eps0')*params_gain.nTHz*Constants('c')/8*(abs((I0_x)*normalization_factor).^2);
        I0_of_x(x_idx) = I0_x*Width*Height*1e3; %mW;
        FWHM_of_x(x_idx) = FWHM_x;
        %
        %     if higher_harmonic
        %
        %         single_rt_data = envelope_x(1:iter_per_rt);
        %         [pks,locs,w,heights] = findpeaks(abs(single_rt_data)/max(abs(single_rt_data)),'MinPeakHeight',0.5);
        %         distance_between_pulses = diff(locs)*dt;
        %
        %     end
        
        
        
    end
 
    ax = plotyy(x(1:skip:end),FWHM_of_x,x(1:skip:end),I0_of_x)
    set(ax(1).YLabel,'String','FWHM (ps)');
    set(ax(2).YLabel,'String','Pulse power (mW)');
    set(ax(1).XLabel,'String','Propagation distance (mm)');
    
    
    ax(1).Children.Color = [088 088  090]/ 255;
    ax(1).Children.LineWidth = 2;
    ax(1).YColor = [088 088  090]/ 255;
    
    ax(2).Children.Color = [227  114  34]/255;
    ax(2).Children.LineWidth = 2;
    ax(2).YColor = [227  114  34]/255;
    
    
    skip_t = 10
    skip_x = 10
    X = x(1:skip:end)';
    X = X(1:skip_x:end);
    figure;
    T_GAIN = [1:skip_t:iter_per_rt]'*dt;
    [X,T_GAIN] = meshgrid(X,T_GAIN);
    surf_data = pulse_data_vs_x(:,1:iter_per_rt)';
    Z = (abs(surf_data(1:skip_t:end,1:skip_x:end)).^2)/max(max(abs(surf_data).^2));
    surf(T_GAIN,X,Z)
    
    figure;
    norm_ = max(abs(F).^2);
    ax = plotyy(x,abs(F).^2/norm_,x,ruu-rll);
    set(ax(1).YLabel,'String','Normalized intensity');
    set(ax(2).YLabel,'String','Population inversion');
    set(ax(1).XLabel,'String','Propagation distance (mm)');
    
    
    ax(1).Children.Color = [088 088  090]/ 255;
    ax(1).Children.LineWidth = 2;
    ax(1).YColor = [088 088  090]/ 255;
    
    ax(2).Children.Color = [162  173  0]/255;
    ax(2).Children.LineWidth = 2;
    ax(2).YColor = [162  173  0]/255;
end
%%
dummy = envelope(intrv1);
NFFT = 2^nextpow2(4*length(dummy));
Y = ifft(dummy.*hanning(length(dummy)),NFFT);
freq = [0:NFFT/2-1, -NFFT/2:-1]*1/dt/NFFT;
[fpks,flocs]  = findpeaks(abs(Y(1:NFFT/2-1)),freq(1:NFFT/2-1),'MinPeakDistance',0.7*1/T_RT); 
frt = flocs(2) - flocs(1); 
TR_ = 1/frt; 

figH = figure;


plot(fftshift(freq), fftshift(abs(Y).^2)/max(abs(Y).^2));
ax = figH.Children;
ax.YLabel.String = 'Optical spectrum (a.u.)';
ax.XLabel.String = 'Envelope frequency (THz)';
ax.XLim = [-.3,.3];
ax.YLim = [0,1.2];


ax.Children.Color = [088 088  090]/ 255;
ax.Children.LineWidth = 2;
ax.YColor = [088 088  090]/ 255;


%% plot the inversion data at a particular position and compare with the analytical solution I had 


%% GAIN 
x_idx_GAIN = 15;
envelope_at_GAIN_X = pulse_data_vs_x(x_idx_GAIN,:);
E0_sqrt_at_GAIN_X = abs(envelope_at_GAIN_X).^2;
E0_sqrt_max_GAIN = max(E0_sqrt_at_GAIN_X);
R_GAIN = gain_model.T2*E0_sqrt_max_GAIN; % E IS NORMALIZED TO mu/hbar E not to mu/2*hbar E!!!! 


d_eq_GAIN = p*d_th;
T1_GAIN = gain_model.T1; 


inversion_at_x_idx_GAIN = inversion_data_gain_vs_x(x_idx_GAIN,:);
times_at_GAIN_X = linspace(0,length(inversion_at_x_idx_GAIN),length(inversion_at_x_idx_GAIN))*dt;


% TR_ = 30.2730; % measured round trip time;

% middle of the pulse; 
% t0_GAIN = 55.4304;  %  for T1g = 20 (ring cavity/noolap) and p=1.3
% t0_GAIN = 60.5760; % for T1g = 20 (ring cavity/noolap) and p=1.2
[pks,locs ] = findpeaks(E0_sqrt_at_GAIN_X/E0_sqrt_max_GAIN,'MinPeakHeight',.6); 
if length(locs)>=2
    t0_GAIN = times_at_GAIN_X(locs(3)); 
else
   t0_GAIN = times_at_GAIN_X(locs(1)); 
end



tau_p_GAIN = FWHM_of_x(x_idx_GAIN)/1.76;
% tau_p = 2.2625;
T_GAIN = tau_p_GAIN*1.76*2; 
h = 2*tau_p_GAIN/T_GAIN;

t_a = t0_GAIN-TR_/2; t_b = t_a + TR_;
t_l = t0_GAIN-T_GAIN/2; t_r = t0_GAIN+T_GAIN/2; 

dT_GAIN = TR_-T_GAIN;

gamma_1_GAIN = 1/T1_GAIN; 
gamma_2_GAIN = h*R_GAIN+gamma_1_GAIN;


dt_small  = 30;
t_idx = times_at_GAIN_X >= t_a-dt_small & times_at_GAIN_X <= t_b+dt_small; 
t_vec_new_GAIN = times_at_GAIN_X(t_idx); 
I_vec_new_GAIN = E0_sqrt_at_GAIN_X(t_idx); 
inv_vec_new_GAIN = inversion_at_x_idx_GAIN(t_idx);


Delta_l_GAIN = d_eq_GAIN*( gamma_1_GAIN/gamma_2_GAIN * (exp(-gamma_2_GAIN*T_GAIN) - 1) - (exp(gamma_1_GAIN*dT_GAIN)-1) )/(exp(-gamma_2_GAIN*T_GAIN)-exp(gamma_1_GAIN*dT_GAIN))
Delta_r_GAIN = (Delta_l_GAIN-d_eq_GAIN*gamma_1_GAIN/gamma_2_GAIN)*exp(-gamma_2_GAIN*T_GAIN)+d_eq_GAIN*gamma_1_GAIN/gamma_2_GAIN
Delta_a_GAIN = (Delta_l_GAIN-d_eq_GAIN)*exp(gamma_1_GAIN*dT_GAIN/2)+d_eq_GAIN


boxcar_pulse_GAIN = h*E0_sqrt_max_GAIN*heaviside(t_vec_new_GAIN - t_l).*heaviside(t_r-t_vec_new_GAIN);




inv_analytical_GAIN = 0*inv_vec_new_GAIN;
d1_idx = t_vec_new_GAIN < t_l;
inv_analytical_GAIN(d1_idx)= (Delta_a_GAIN-d_eq_GAIN)*exp(-gamma_1_GAIN*(t_vec_new_GAIN(d1_idx)-t_a))+d_eq_GAIN; 

d2_idx = t_vec_new_GAIN >= t_l & t_vec_new_GAIN <= t_r;
inv_analytical_GAIN(d2_idx)= (Delta_l_GAIN-d_eq_GAIN*gamma_1_GAIN/gamma_2_GAIN)*exp(-gamma_2_GAIN*(t_vec_new_GAIN(d2_idx)-t_l)) ...
    + d_eq_GAIN*gamma_1_GAIN/gamma_2_GAIN; 

d3_idx = t_vec_new_GAIN > t_r;
inv_analytical_GAIN(d3_idx)= (Delta_r_GAIN-d_eq_GAIN)*exp(-gamma_1_GAIN*(t_vec_new_GAIN(d3_idx)-t_r))+d_eq_GAIN; 

figure;
dTH_vector = d_th*(1+0*t_vec_new_GAIN)
plotyy(t_vec_new_GAIN,[I_vec_new_GAIN;boxcar_pulse_GAIN], ...
    t_vec_new_GAIN,[inv_vec_new_GAIN;inv_analytical_GAIN;dTH_vector]); 


trapz(t_vec_new_GAIN,inv_vec_new_GAIN)
trapz(t_vec_new_GAIN,inv_analytical_GAIN)

%%%% the following is from the analytical expression for the integral! 
int_1 = (Delta_a_GAIN-d_eq_GAIN)/gamma_1_GAIN* ... 
    (1-exp(-gamma_1_GAIN*dT_GAIN/2))+dT_GAIN*d_eq_GAIN/2;

int_2 = 1/gamma_2_GAIN*(Delta_l_GAIN-d_eq_GAIN*gamma_1_GAIN/gamma_2_GAIN)* ...
    (1-exp(-gamma_2_GAIN*T_GAIN))+T_GAIN*d_eq_GAIN*gamma_1_GAIN/gamma_2_GAIN;

int_3 = (Delta_r_GAIN-d_eq_GAIN)/gamma_1_GAIN* ... 
    (1-exp(-gamma_1_GAIN*dT_GAIN/2))+dT_GAIN*d_eq_GAIN/2;





%% ABSORBER 
offset = size(inversion_data_gain_vs_x,1);
x_idx_ABS = 10;

envelope_at_ABS_X = pulse_data_vs_x(x_idx_ABS+offset,:);
E0_sqrt_at_ABS_X = abs(envelope_at_ABS_X).^2;
E0_sqrt_max_ABS = max(E0_sqrt_at_ABS_X);
R_ABS = (abs_model.zUL/abs_model.zNORM)^2*abs_model.T2*E0_sqrt_max_ABS; % E IS NORMALIZED TO mu/hbar E not to mu/2*hbar E!!!! 


d_eq_ABS = -1;
T1_ABS = abs_model.T1; 


TR_ = 30.2730; % measured round trip time; --> automate


inversion_at_x_idx_ABS = inversion_data_abs_vs_x(x_idx_ABS,:);
% mid of the pulse; --> automate
% t0_ABS = 49.2136;  %  for T1g = 20 (ring cavity/noolap) and p=1.2
% t0_ABS = 54.4793;  %  for T1g = 20 (ring cavity/noolap) and p=1.3
times_at_ABS_X = linspace(0,length(inversion_at_x_idx_ABS),length(inversion_at_x_idx_ABS))*dt;
[pks,locs ] = findpeaks(E0_sqrt_at_ABS_X/E0_sqrt_max_ABS,'MinPeakHeight',.6); 
if length(locs)>=2
    t0_ABS = times_at_ABS_X(locs(2)); 
else
   t0_ABS = times_at_ABS_X(locs(1)); 
end



tau_p = FWHM_of_x(x_idx_ABS+offset)/1.76;
T_ABS = tau_p/.50; 
h = 2*tau_p/T_ABS;

t_a = t0_ABS-TR_/2; t_b = t_a + TR_;
t_l = t0_ABS-T_ABS/2; t_r = t0_ABS+T_ABS/2; 

dT_ABS = TR_-T_ABS;

gamma_1_ABS = 1/T1_ABS; 
gamma_2_ABS = h*R_ABS+gamma_1_ABS;


 
t_idx = times_at_ABS_X >= t_a & times_at_ABS_X <= t_b; 
t_vec_new_ABS = times_at_ABS_X(t_idx); 
I_vec_new_ABS = E0_sqrt_at_ABS_X(t_idx); 
inv_vec_new_ABS = inversion_at_x_idx_ABS(t_idx);


Delta_l_ABS = d_eq_ABS*( gamma_1_ABS/gamma_2_ABS * (exp(-gamma_2_ABS*T_ABS) - 1) - (exp(gamma_1_ABS*dT_ABS)-1) )/(exp(-gamma_2_ABS*T_ABS)-exp(gamma_1_ABS*dT_ABS))
Delta_r_ABS = (Delta_l_ABS-d_eq_ABS*gamma_1_ABS/gamma_2_ABS)*exp(-gamma_2_ABS*T_ABS)+d_eq_ABS*gamma_1_ABS/gamma_2_ABS
Delta_a_ABS = (Delta_l_ABS-d_eq_ABS)*exp(gamma_1_ABS*dT_ABS/2)+d_eq_ABS


boxcar_pulse_ABS = h*E0_sqrt_max_ABS*heaviside(t_vec_new_ABS - t_l).*heaviside(t_r-t_vec_new_ABS);




inv_analytical_ABS = 0*inv_vec_new_ABS;
d1_idx = t_vec_new_ABS < t_l;
inv_analytical_ABS(d1_idx)= (Delta_a_ABS-d_eq_ABS)*exp(-gamma_1_ABS*(t_vec_new_ABS(d1_idx)-t_a))+d_eq_ABS; 

d2_idx = t_vec_new_ABS >= t_l & t_vec_new_ABS <= t_r;
inv_analytical_ABS(d2_idx)= (Delta_l_ABS-d_eq_ABS*gamma_1_ABS/gamma_2_ABS)*exp(-gamma_2_ABS*(t_vec_new_ABS(d2_idx)-t_l)) ...
    + d_eq_ABS*gamma_1_ABS/gamma_2_ABS; 

d3_idx = t_vec_new_ABS > t_r;
inv_analytical_ABS(d3_idx)= (Delta_r_ABS-d_eq_ABS)*exp(-gamma_1_ABS*(t_vec_new_ABS(d3_idx)-t_r))+d_eq_ABS; 

figure;
plotyy(t_vec_new_ABS,[I_vec_new_ABS;boxcar_pulse_ABS],t_vec_new_ABS,[inv_vec_new_ABS;inv_analytical_ABS]); 


trapz(t_vec_new_ABS,inv_vec_new_ABS)
trapz(t_vec_new_ABS,inv_analytical_ABS)



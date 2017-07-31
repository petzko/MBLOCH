close all

% dt = dat.dt;
envelope = record_F_a;


iter_per_rt = round(T_RT/dt/2);

lbdry = 500*iter_per_rt;
rbdry = 508*iter_per_rt;
intrv1 = [lbdry:rbdry];
% this is the simulation data intensity normalized to 1 for max intensity
data1_unnormalized = abs(envelope(intrv1)).^2;
data1 = data1_unnormalized/max(data1_unnormalized);
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
p = (gain_model.rho_u_DC_eq-gain_model.rho_l_DC_eq)/d_th

gth = d_th*sigma_a*Ncarriers_a*params_abs.L*1e-3
g0 = p*gth


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

q0 =  p*d_th*sigma_g*Ncarriers_g*params_gain.L*1e-3;
aLg = 2*params_gain.linear_loss*100*params_gain.L*1e-3;

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
        [pks,locs,w,heights] = findpeaks(abs(envelope_x)/max(abs(envelope_x)),'MinPeakHeight',0.01);
        
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
    T = [1:skip_t:iter_per_rt]'*dt;
    [X,T] = meshgrid(X,T);
    surf_data = pulse_data_vs_x(:,1:iter_per_rt)';
    Z = (abs(surf_data(1:skip_t:end,1:skip_x:end)).^2)/max(max(abs(surf_data).^2));
    surf(T,X,Z)
    
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
figure; 
plot(x,abs(F).^2,'-b',x,abs(B).^2);

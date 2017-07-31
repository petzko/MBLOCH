close all

Width = 40*1E-6; % 20um 
Height = 10e-6; 

n = 3.6;
tch=1e-12;
zUL = sim_settings.zULa;
dt = dat.dt;
envelope = record_U_a;


normalization_factor = 1E12*Constants('hbar')/(Constants('q0')*zUL*1E-9);
It = Constants('eps0')*n*Constants('c')/8*(abs((envelope)*normalization_factor).^2); % in watts
pulse_power = It*Width*Height*1e3;


Ix_ = Constants('eps0')*n*Constants('c')/8*(abs((dat.U)*normalization_factor).^2);
Px_ = Ix_*Width*Height*1e3;

time =linspace(0, dt*length(It),length(It))*tch;


 
tms = [0:length(record_U_g)-1]*dat.dt;

lbdry = 280*dat.T_R
rbdry = 296*dat.T_R


win = @hanning; 
NFFT = 2^nextpow2(length(envelope)); 
interval = tms > 100*dat.T_R & tms < 250*dat.T_R; 
Y = ifft(envelope.*win(length(envelope)),NFFT);
f = 1/NFFT*[0:NFFT/2-1, -NFFT/2:-1]*1/dat.dt;

figsize = [0,0,0.2,0.2]
figure('units','normalized','position',figsize)
ax = plotyy(tms/dat.T_R,pulse_power,tms/dat.T_R,record_r22g-record_r11g);

set(ax(2).Children(1),'Color','r');
set(ax(1).Children(1),'Color','b');
set(ax,{'ycolor'},{'b';'r'})
ylabel(ax(1), 'P(t) (mW)');
ylabel(ax(2), 'pop. inversion');
xlabel('time (T_{rt})'); 
xlim(ax(1),[lbdry,rbdry]/dat.T_R)
xlim(ax(2),[lbdry,rbdry]/dat.T_R)

idx_t = tms > lbdry & tms < rbdry;
pulses = abs(envelope(idx_t)).^2/max(abs(envelope(idx_t)).^2);
times_t = tms(idx_t);


diff_T = times_t(2)-times_t(1);
idx_s = pulses > 0.49 & pulses < 0.51; 
pulsewidths = times_t(idx_s);
widths = diff(pulsewidths);
width = unique(widths(widths > 1 & widths < 1/2*dat.T_R));
display(['pulse width : ' num2str(width)])


figure('units','normalized','position',figsize)
subplot(1,2,2); 
ax = plotyy(dat.x,Px_,dat.x,r22-r11);
set(ax(2).Children(1),'Color','r');
set(ax(1).Children(1),'Color','b');

set(ax,{'ycolor'},{'b';'r'})
xlabel('Length (mm)'); 
ylabel(ax(1), 'Power (mW)');
ylabel(ax(2), 'Inversion');

subplot(1,2,1); 
ax = plot(f,abs(Y).^2,'b');
xlabel('Freq. (THz)'); ylabel('|E(\omega)|^2 (a.u.)');
xlim([-.4,.4])


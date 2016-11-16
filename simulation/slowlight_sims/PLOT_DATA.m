close all
tms = [0:length(record_U_g)-1]*dat.dt;

lbdry = 190*dat.T_R
rbdry = 195*dat.T_R


win = @hanning; 
NFFT = 2^nextpow2(length(record_U_g)); 
interval = tms > 100*dat.T_R & tms < 250*dat.T_R; 
Y = ifft(record_U_g.*win(length(record_U_g)),NFFT);
f = 1/NFFT*[0:NFFT/2-1, -NFFT/2:-1]*1/dat.dt;


subplot(2,2,1); 
ax = plotyy(tms,abs(record_U_g).^2,tms,record_r22g-record_r11g);
title('Gain'); xlabel('time (ps)'); 
set(ax(2).Children(1),'Color','r');
set(ax(1).Children(1),'Color','b');
set(ax,{'ycolor'},{'b';'r'})
ylabel(ax(1), 'I(t)');
ylabel(ax(2), 'pop. inversion');
xlim(ax(1),[lbdry,rbdry])
xlim(ax(2),[lbdry,rbdry])

subplot(2,2,2); 
ax = plotyy(tms,abs(record_U_a).^2,tms,record_r22a-record_r11a);
title('Absorber'); xlabel('time (ps)'); 
set(ax(2).Children(1),'Color','r');
set(ax(1).Children(1),'Color','b');
set(ax,{'ycolor'},{'b';'r'})
ylabel(ax(1), 'I(t)');
ylabel(ax(2), 'pop. inversion');
xlim(ax(1),[lbdry,rbdry]);
xlim(ax(2),[lbdry,rbdry]);

subplot(2,2,[3,4]); 
ax = plot(f,abs(Y).^2,'b');
xlabel('freq. THz'); ylabel('|E(\omega)|^2');
xlim([-1.5,1.5])



figure; 
ax = plotyy(dat.x,abs(dat.U).^2,dat.x,r22-r11);
set(ax(2).Children(1),'Color','r');
set(ax(1).Children(1),'Color','b');

set(ax,{'ycolor'},{'b';'r'})
xlabel('LENGTH (mm)'); 
ylabel(ax(1), 'I(t)');
ylabel(ax(2), 'pop. inversion');
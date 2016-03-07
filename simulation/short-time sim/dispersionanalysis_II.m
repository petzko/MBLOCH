close all;

trip_1 = 1; Pt_1 =2; signal_1 = squeeze(U_t(Pt_1,trip_1,:)); 
trip_2 = 1 ; Pt_2 =5; signal_2 = squeeze(U_t(Pt_2,trip_2,:));

dL = Ltot*(trip_2-trip_1) + x(checkpoints(Pt_2))-x(checkpoints(Pt_1)); 

times = linspace(0,1,length(signal_1))'*dt*recordIter*length(signal_1); 

win_1 = (hanning(length(signal_1))); 
signal_1 = signal_1.*win_1; 
win_2 = (hanning(length(signal_2))); 
signal_2 = signal_2.*win_2; 

NFFT = 2^nextpow2(32*length(signal_1));
% NFFT = length(signal_1);
Y1 = ifft(signal_1,NFFT); Y1 = fftshift(Y1); 
Y2 = ifft(signal_2,NFFT); Y2 = fftshift(Y2);

f0 = E0/2/pi; fs = 1/dt; 
f1 = linspace(-1/2,1/2,length(Y1))'*fs+f0; df = f1(2)- f1(1);  dw = 2*pi*df;

lim_idx = and(f1 > f0-.7, f1<=f0+.7);
f1 = f1(lim_idx);w1 = 2*pi*f1;
Y1 = Y1(lim_idx); Y2 = Y2(lim_idx); 
Psi_01 = (unwrap((angle(Y2)-angle(Y1))))+w1/c*dL; 
Psi_02 = unwrap((angle(Y2./Y1))); 
Psi_03 = angle(Y2./Y1);

gn = log((abs(Y2)./abs(Y1)));
w1 = f1*2*pi;

n_R = c_0.*Psi_01./(f_1*2*pi*dL)+n*f0./(f_1); 
n_Im = c_0.*gn./(f1*2*pi*dL);
%wave number group velocity GVD and GDD. 
k_w = n_R*2*pi.*f1/c_0;

GVD = diff(diff(k_w))/dw^2; 

short_idx = and(f1 > f0-10*df,f1 <f0+10*df); 
k_w_short = detrend(k_w(short_idx),'linear');  f_1_short=  f1(short_idx)-f0;


figure;
subplot(1,3,1); 
lim = [3.2, 4.6];
plot(f1,n_R); ylabel('\Re\{n\}'); xlabel('Freq. [THz]');xlim(lim); 
subplot(1,3,2);    
plot(f1,n_Im); ylabel('\Im\{n\}'); xlabel('Freq. [THz]');xlim(lim); 
subplot(1,3,3);
plot(f1,n_Im.*f1./c_0*10); ylabel('gain (1/cm)'); xlabel('Freq. [THz]');xlim(lim); 


figure;
subplot(1,2,1); 
plot(f1,k_w/2/pi,f1,f1*n/c_0); ylabel('wave number (1/mm)'); xlabel('Freq. [THz]');xlim(lim); 
subplot(1,2,2); 
plot(f1(2:end-1),GVD); ylabel('Group Velocity Dispersion (ps^2/mm)'); xlabel('Freq. [THz]');xlim(lim); 

% figure; 
% w_1_short = 2*pi*f_1_short;
% 
% k_j = polyfit(w_1_short,k_w_short,3);
% k_3 = 6*k_j(1);
% k_2 = 2*k_j(2);
% v_g = k_j(3) ;
% k_0 = k_j(4);
% % remove residual linear wavenumber!
% k_w_short_new = k_w_short - v_g*w_1_short-k_0;
% plot(w_1_short,k_w_short_new);
% k_j = polyfit(w_1_short,k_w_short_new,3); ylabel('detrended wave number (1/mm)'); xlabel('Freq. [THz]'); 
% k_3 = 6*k_j(1)
% k_2 = 2*k_j(2)
% v_g = k_j(3) 
% k_0 = k_j(4) 
% 
% plot(f_1_short,k_w_short,f_1_short,k_w_short_new); ylabel('detrended wave number (1/mm)'); xlabel('Freq. [THz]'); 

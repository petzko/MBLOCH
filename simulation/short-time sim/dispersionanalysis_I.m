close all;

trip_1 = 1; Pt_1 =1; signal_1 = squeeze(U_t(Pt_1,trip_1,:)); 
trip_2 = 3; Pt_2 =1; signal_2 = squeeze(U_t(Pt_2,trip_2,:));

dL = Ltot*(trip_2-trip_1) + x(checkpoints(Pt_2))-x(checkpoints(Pt_1)); 
times = linspace(0,1,length(signal_1))'*dt*recordIter*length(signal_1); 

options = {'hanning','truncate'};
win_1 = [hanning(length(signal_1))]; 
signal_1 = signal_1.*win_1; 
win_2 = (hanning(length(signal_2))); 
signal_2 = signal_2.*win_1; 


subplot(1,2,1); 
plot(times,real(signal_1.*exp(-1i*times*E0)).*win_1);
subplot(1,2,2); 
plot(times,real(signal_2.*exp(-1i*times*E0)).*win_2);


% Y_1 = fft(real((signal_1).*exp(-1i*times*E0)));
% Y_2 = fft(real((signal_2).*exp(-1i*times*E0))); 
% f_1 = 1/dt*linspace(0,1,length(Y_1)); f_2 = f_1;
b0 = E0*n/c_0;
[f_1, Y_1] = FourierTransf( real(signal_1.*exp(-1i*times*E0)),dt,0,1,1,options);
[f_2, Y_2] = FourierTransf( real(signal_2.*exp(-1i*times*E0)),dt,0,1,1,options);
df = f_1(2)-f_1(1);
% 
% figure; plotyy(f_1,abs(Y_1),f_2,abs(Y_2));


A_w = (Y_2./Y_1);  Psi = (angle(A_w)); %in rad  
f0 = E0/2/pi; 

n_Re = -c_0/dL./f_1'/2/pi.*Psi + n*f0./f_1';  
n_Im = c_0/dL./f_1'/2/pi.*log((abs(Y_2)./abs(Y_1)));



GVD_1 = 2/c_0*diff(n_Re)/df;
GVD_2 = f_1(2:end-1)'./c_0/df/df.*diff(diff(n_Re));
GVD_m2 = GVD_1(2:end)+GVD_2; 

lim = [f0-0.5,f0+0.5];

beta = f_1'.*n_Re./c_0;
b_idx = and(f_1' > E0/2/pi-5*df,f_1' <E0/2/pi+5*df); beta_short = detrend(beta(b_idx));  f_1_short=  f_1(b_idx)-E0/2/pi; 
figure; plot(f_1_short,(beta_short)); xlabel(['Freq. relative to f_0 = ' num2str(E0/2/pi) ' [THz]']); ylabel('\beta(f) 1/mm');




figure;
v_g = df./(diff(beta)); beta_2 = diff(diff(beta))/df^2;  beta_3 = diff(diff(diff(beta)))/df^3;
subplot(1,2,1); plot(f_1,n_Re);xlim(lim); xlabel('Freq. [THz]'); ylabel('\Re(n(f))');  legend(num2str(ampl));
subplot(1,2,2); plot(f_1,n_Im);xlim(lim); xlabel('Freq. [THz]'); ylabel('\Im(n(f))');  legend(num2str(ampl));

% 
figure; 
subplot(1,3,1); plot(f_1(2:end),v_g); xlim(lim); xlabel('Freq. [THz]'); ylabel('v_{gr}(f)');
subplot(1,3,2); plot(f_1(1:end-2),[beta_2,GVD_m2]); xlim(lim); xlabel('Freq. [THz]'); ylabel('\beta_2 ps^2/mm');
subplot(1,3,3); plot(f_1(1:end-3),beta_3); xlim(lim); xlabel('Freq. [THz]'); ylabel('\beta_3 ps^2/mm');
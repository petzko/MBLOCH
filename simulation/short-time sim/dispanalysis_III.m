close all; 
f0= E0/2/pi; 

%%%% collect the signals
trip_1 = 1; Pt_1 =2; signal1 = squeeze(U_t(Pt_1,trip_1,:)); signal1 = padarray(signal1,1E3,0,'pre');
 tms1 = linspace(0,dt*length(signal1),length(signal1))'; 
trip_2 = 1 ; Pt_2 =4; signal2 = squeeze(U_t(Pt_2,trip_2,:));signal2 = padarray(signal2,1E3,0,'pre');
tms2 = linspace(0,dt*length(signal2),length(signal2))'; 



%calculate the propagated distance
dL = Ltot*(trip_2-trip_1) + x(checkpoints(Pt_2))-x(checkpoints(Pt_1)); 

%calculate the central freq's propagation time
avg_delay = dL/c;

%plot the envelopes
figure;  subplot(1,2,1);  plot(tms1,real(signal1.*exp(-1i*E0*tms1))); xlabel('time (ps)'); ylabel('e_in(t)'); 

subplot(1,2,2);  plot(tms2,real(signal2.*exp(-1i*E0*tms2))); xlabel('time (ps)') ; ylabel('e_out(t)');

win = hanning(length(signal1))*0+1; 
NFFT = 2^nextpow2(16*length(signal1))*0+length(signal1);
% fft transform the envelopes and center them. 
Y1 = ifft(signal1.*win,NFFT); Y1 = fftshift(Y1);
Y2 = ifft(signal2.*win,NFFT); Y2 = fftshift(Y2); 

%take the phase
Psi = unwrap(angle(Y2./Y1));
%frequency.. 
f = linspace(-1/2,1/2,length(Y1))'*1/dt; df = f(2)-f(1); 
% figure; subplot(1,2,1); plot(f,abs(Y1));  subplot(1,2,2); plot(f,abs(Y2)); 

%bandwidth of the measurement
BW = .8; %THz
lim_idx= and(f>-BW/2,f<=BW/2);
%extract the phase within the desired bandwidth. 
f_ = f(lim_idx);w_ =2*pi*f_; 
Psi_ = (Psi(lim_idx));

% remove the constant part of the acquired phase
zer_idx = and(f_>=(-df-2*eps),f_<=(df+2*eps));
psi_0 = mean(Psi_(zer_idx));
Psi_ = Psi_ - psi_0; 

%estimate the gain:
gn = log(abs(Y2)./abs(Y1))./dL;
gn = gn(lim_idx);

% get the real refractive index:
n_Re = c_0/dL*Psi_./(2*pi*(f_ +f0)) + n*f0./(f_+f0);
n_Im =c_0./(2*pi*(f_+f0)).*gn; 

beta_w = n_Re.*(w_+E0)/c_0; 
vg = dw./diff(beta_w);
GVD = diff(diff(beta_w))/dw^2; 

figure;
subplot(1,4,1); plot(f_+f0,n_Re); ylabel('n_{Re}'); xlabel('Freq. (THz)');
subplot(1,4,2); plot(f_+f0,n_Im); ylabel('n_{Im}'); xlabel('Freq. (THz)');
subplot(1,4,3); plot(f_(1:end-1) +f0,vg/c); ylabel('v_{gr}/v_ph )'); xlabel('Freq. (THz)');
subplot(1,4,4); plot(f_(2:end-1) +f0,GVD); ylabel('GVD (ps^2/mm)'); xlabel('Freq. (THz)');


betas = polyfit(w_,beta_w,3); 
b1 = betas(3);
b0 = betas(4);
beta_w = beta_w - betas(4)-betas(3).*w_;
betas = polyfit(w_,beta_w,3); 

b3 = betas(1)*6
b2 = betas(2)*2


figure; 
plot(f_,beta_w,f_,betas(1).*w_.^3+betas(2).*w_.^2+betas(4).*w_+betas(4));
xlabel('Freq. (THz)'); ylabel('\beta(\omega) and \beta_{approx}(\omega)'); 
tau_3 = dL*(b1+b2*w_+b3/2*w_.^2);
tau_2 = dL*(b1+b2*w_);

figure; plot(f_,tau_2,f_,tau_3); 




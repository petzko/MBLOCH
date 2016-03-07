close all;
f0= E0/2/pi;
plotON = true;
convention = 'physics';
if(strcmp(convention,'physics'))
    transform = @ifft;
else
    transform = @fft;
end

%%%% collect the signals
trip_1 = 1; Pt_1 =2; 
signal1 = squeeze(U_t(Pt_1,:))';
% signal1=signal1(1:length(signal1)/2);
trip_2 = 1 ; Pt_2 =3;
signal2 = squeeze(U_t(Pt_2,:))';
signal1 = padarray(signal1,1E5,'both');
signal2 = padarray(signal2,1E5,'both');


tms1 =[0:length(signal1)-1].'*dt;
tms2 =[0:length(signal2)-1].'*dt;

%calculate the propagated distance
dL = Ltot*(trip_2-trip_1) + x(checkpoints(Pt_2))-x(checkpoints(Pt_1));

%calculate the central freq's propagation time
avg_delay = dL/c;

% calculate the normalization constant for the fourier transform!
win = hanning(length(signal1));
NFFT = 1E6;
nconst = normalizefourier(NFFT,win,transform);

%plot the envelopes
if plotON
%     figure;
    dfigure('DName','Pre and post-interaction pulses');
    normconst1 = max(abs(signal1)); 
    normconst2 = max(abs(signal2)); 
    subplot(2,2,1);  
    plot(tms1,real(signal1.*exp(-1i*E0*tms1))/normconst1); xlabel('time (ps)'); ylabel('E_{in}(t)'); set(gca, 'box', 'off'); 
    subplot(2,2,2);  
    ax = plot(tms2,real(signal2.*exp(-1i*E0*tms2)/normconst1)); xlabel('time (ps)'); ylabel('E_{out}(t)');
    set(gca, 'box', 'off')
end


% Fourier transform the envelopes and center them.
Y1 = transform(signal1.*win,NFFT)/nconst; Y1 = fftshift(Y1);
Y2 = transform(signal2.*win,NFFT)/nconst; Y2 = fftshift(Y2);

% Y1 = transform(signal1)/nconst; Y1 = fftshift(Y1);
% Y2 = transform(signal2)/nconst; Y2 = fftshift(Y2);

%take the phase
Psi = unwrap(angle(Y2)-angle(Y1));
%calculate the frequency..
fs = 1/dt; f = fftshift([[0:length(Y1)/2-1],[-length(Y1)/2:-1]].')*fs/length(Y1); df = f(2)-f(1); %in THz

%allowed bandwidth of the measurement
BW = 1.0; %THz
lim_idx= and(f>-BW/2,f<=BW/2);
%extract the phase within the desired bandwidth.
f_ = f(lim_idx); w_ =2*pi*f_; dw = w_(2)-w_(1);
Psi_ = (Psi(lim_idx));

% remove any constant part of the acquired phase
zer_idx = and(f_>=(-df-2*eps),f_<=(df+2*eps));
psi_0 = mean(Psi_(zer_idx)); Psi_ = Psi_ - psi_0;

%estimate the gain:
gn = log(abs(Y2)./abs(Y1))./dL;
gn = gn(lim_idx);

% get the real refractive index:
n_Re = c_0/dL*Psi_./(2*pi*(f_ +f0)) + n*f0./(f_+f0);
n_Im = -c_0./(2*pi*(f_+f0)).*gn;

%get the wave number, the group velocity and the GVD
beta_w_tot = n_Re.*(w_+E0)/c_0;
vg = dw./diff(beta_w_tot);
GVD = diff(diff(beta_w_tot))/dw^2;
Nw = length(beta_w_tot);Nw = Nw-mod(Nw,2);
beta_w = beta_w_tot-beta_w_tot(Nw/2)-(beta_w_tot(Nw/2+1)-beta_w_tot(Nw/2))*w_/dw; 
%plot the results so far

BW2 = .7; %THz
lim_idx2= and(f_>-BW2/2,f_<=BW2/2);

betas = polyfit(w_(lim_idx2),beta_w(lim_idx2),3);
beta_w = beta_w - betas(4)-betas(3).*w_;

betas = polyfit(w_(lim_idx2),beta_w(lim_idx2),3);
beta_poly = betas(1).*w_.^3+betas(2).*w_.^2+betas(4).*w_+betas(4);
%calculate the second order and third order dispersion coefficients.
b3 = betas(1)*6
b2 = betas(2)*2
b1 = betas(3)
b0 = betas(4)

if (plotON)
    subplot(2,2,[3 4]);
    plot(f_+f0,beta_w,f_+f0,beta_poly);
    xlabel('Freq. (THz)'); ylabel('\beta(\omega) and \beta\prime(\omega)'); 
    
end
% close all;
hbar = Constants('hbar',{'time',params_abs.tch})/Constants('q0');
E0 = params_gain.E0/hbar;

f0= E0/2/pi;
plotON = true;
convention = 'physics';
if(strcmp(convention,'physics'))
    transform = @ifft;
else
    transform = @fft;
end

analyze_gain = 1;
analyze_abs = 1;
analyze_all = analyze_gain*analyze_abs;

if analyze_gain && ~analyze_abs
    %%%% collect the signals
    signal1 = record_U_G_IN;
    signal2 = record_U_G_OUT;
    dL = params_gain.L;
    
elseif analyze_abs && ~analyze_gain
    %%%% collect the signals
    signal1 = record_U_A_IN;
    signal2 = record_U_A_OUT;
    dL = params_abs.L;
else
    signal1 = record_U_G_IN;
    signal2 = record_U_A_OUT;
    dL = params_abs.L+params_gain.L;
end



signal1 = padarray(signal1,1E3,'both');
signal2 = padarray(signal2,1E3,'both');

tms1 =[0:length(signal1)-1].'*dt;
tms2 =[0:length(signal2)-1].'*dt;

%calculate the propagated distance

%calculate the central freq's propagation time

% the effective index 
neff = (params_gain.nTHz*params_gain.L/Ltot+params_abs.nTHz*params_abs.L/Ltot);
c_0 = c*neff;
avg_delay = dL/c;

% calculate the normalization constant for the fourier transform!
win = hanning(length(signal1));
NFFT = 1E6;
%
% win = ones(length(signal1),1);
% NFFT = length(signal1);

nconst = normalizefourier(NFFT,win,transform);

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
BW = 3.5; %THz
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
n_Re = c_0/dL*Psi_./(2*pi*(f_ +f0)) + neff*f0./(f_+f0);
n_Im = -c_0./(2*pi*(f_+f0)).*gn;

%get the wave number, the group velocity and the GVD
beta_w = n_Re.*(w_+E0)/c_0; vg = dw./diff(beta_w); 
GVD = diff(diff(beta_w))/dw^2; 
GDD = GVD*dL;

%plot the results so far
if (plotON)
    
    plotlimits = [-.5,.5]; 
    
    %     subplot(2,3,1); plot(f_+f0,n_Re); ylabel('n_{Re}'); xlabel('Freq. (THz)'); xlim(plotlimits);
    subplot(2,3,[1 2 3]); 
    plot(f_(2:end-1),GDD); hold on;
    ylabel('GDD (\omega) (ps^2)'); xlabel('f-f0 (THz)');
    xlim(plotlimits);
    subplot(2,3,[4 5 6]);
    ax = plotyy(f_,gn*10,f_(1:end-1),vg/c);  
    hold on;
    xlabel('f-f0 (THz)');xlim(plotlimits);
    set(get(ax(1),'Ylabel'),'String','Gain (1/cm)');
    set(get(ax(2),'Ylabel'),'String','v_{gr}/v_{ph}');
    set(ax(1),'xlim',plotlimits);
    set(ax(2),'xlim',plotlimits);
    set(ax(2),'ylim',[.5 , 1.5]);
    
end

% %polynomial fit of the wavenumber
% BW2 = .7; %THz
% lim_idx2= and(f_>-BW2/2,f_<=BW2/2);
% betas = polyfit(w_(lim_idx2),beta_w(lim_idx2),3);
% b1 = betas(3);
% b0 = betas(4);
% %remove the nonlinear part of the wavenumber (i.e. trivial group delay)
% beta_w = beta_w - betas(4)-betas(3).*w_;
% betas = polyfit(w_(lim_idx2),beta_w(lim_idx2),3);
% 
% beta_poly = betas(1).*w_.^3+betas(2).*w_.^2+betas(4).*w_+betas(4);
% %calculate the second order and third order dispersion coefficients.
% b3 = betas(1)*6;
% b2 = betas(2)*2;
% 
% 
% if (plotON)
%     %     dfigure;
%     subplot(2,3,1);
%     plot(f_+f0,beta_w,f_+f0,beta_poly);
%     xlabel('Freq. (THz)'); ylabel('\beta(\omega) and \beta_{fit}(\omega)');xlim([3.2,4.55])
% end
% 
% % plot the group delay with and without 3rd order dispersion.
% tau_3 = dL*(b1+b2*w_+b3/2*w_.^2);
% tau_2 = dL*(b1+b2*w_);
% 
% 

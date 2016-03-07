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
% signal1 = padarray(signal1,1E5,'both');
% signal2 = padarray(signal2,1E5,'both');


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
    [peak1,loc1,width1] = findpeaks(abs(signal1/normconst1),tms1,'MinPeakHeight',0.9);
    [peak2,loc2,width2] = findpeaks(abs(signal2)/normconst2,tms1,'MinPeakHeight',0.9);
    
    
    
    subplot(1,2,1);  
    plot(tms1,real(signal1.*exp(-1i*E0*tms1))/normconst1); xlabel('time (ps)'); ylabel('E_{in}(t)'); set(gca, 'box', 'off'); 
%     xlim([mean(loc1)-3*mean(width1),mean(loc1)+3*mean(width1)]);
    subplot(1,2,2);  
    ax = plot(tms2,real(signal2.*exp(-1i*E0*tms2)/normconst1)); xlabel('time (ps)'); ylabel('E_{out}(t)');
    set(gca, 'box', 'off')
%     xlim([mean(loc2)-3*mean(width2),mean(loc2)+3*mean(width2)]);

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

%plot the results so far
if (plotON)
    plotlimits = [f_(1)+f0,f_(end)+f0]; dfigure;
    subplot(1,5,1); plot(f_+f0,n_Re); ylabel('n_{Re}'); xlabel('Freq. (THz)'); xlim(plotlimits);
    subplot(1,5,2); plot(f_+f0,n_Im); ylabel('n_{Im}'); xlabel('Freq. (THz)');xlim(plotlimits);
    subplot(1,5,3); plot(f_+f0,gn*10); ylabel('gain (1/cm)'); xlabel('Freq. (THz)');%xlim(plotlimits);
    subplot(1,5,4); plot(f_(1:end-1) +f0,vg/c); ylabel('v_{gr}/v_{ph}'); xlabel('Freq. (THz)');xlim(plotlimits);
    subplot(1,5,5); plot(f_(2:end-1) +f0,GVD); ylabel('GVD (ps^2/mm)'); xlabel('Freq. (THz)');xlim(plotlimits);
end
%polynomial fit of the wavenumber
BW2 = .7; %THz
zer_idx = and(f_>=-2*df,f_<=2*df);

slope = mean(diff(beta_w_tot(zer_idx))/dw);
beta_0 = mean(beta_w_tot(zer_idx));
beta_w  = beta_w_tot -slope*w_-beta_0;
% [pks,locs,widths] = findpeaks(gn*10,f_+f0,'MinPeakHeight',0); f1 = locs(1); f2 = locs(2);
[pks,locs,widths] = findpeaks(abs(beta_w),f_+f0,'MinPeakHeight',0); f1 = locs(1); f2 = locs(2);
width = min(widths); 
idx1 = (f_+f0 >=locs(2)-widths(2)/2) & (f_+f0 <=locs(2)+widths(2)/2);
idx2 = (f_+f0 >=locs(5)-widths(5)/2) & (f_+f0 <=locs(5)+widths(5)/2);


poly1 = polyfit(w_(idx1),beta_w(idx1),2);
poly2 = polyfit(w_(idx2),beta_w(idx2),2);

% beta_quad1 = poly1(1)*w_.^3+poly1(2)*w_.^2+poly1(3)*w_+poly1(4);
% beta_quad2 = poly2(1)*w_.^3+poly2(2)*w_.^2+poly2(3)*w_+poly2(4);

beta_quad1 = poly1(1)*w_.^2+poly1(2)*w_+poly1(3);
beta_quad2 = poly2(1)*w_.^2+poly2(2)*w_+poly2(3);

dfigure; plot(f_+f0,beta_w,f_+f0,beta_quad1,f_+f0,beta_quad2); xlabel('Freq. (THz)');
ylabel('Wavenumber variation (1/mm)');
xlim([3.4,4.35]); ylim([-2.5,2.5]);
legnfo = {'\beta(\omega)','\beta^{I}(\omega)','\beta^{II}(\omega)'};
dlegend(legnfo,'Wave number:');

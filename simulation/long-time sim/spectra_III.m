%% clean clear and close stuff
close all; clc; 

convention = 'physics';

%determine the transform function according to the assumed convention 
if(strcmp(convention,'physics'))
    transform = @ifft;
    itransform = @fft; 
else
    transform = @fft;
    itrasform = @ifft;
end

winfunc = @hanning; 
NFFTfunc = @(elem) 2^nextpow2(4*elem); 

load('disp_compensated');

c = 0.299792458/3.6; % central freq phase vel.  
T_R = 2*5/c; %round trip time in a 5mm cavity
f_R = 1/T_R; 
Ntau = 5*round(T_R/Dt);

E0 = 3.88*2*pi; % central frequency (i.e. energy difference between ULL and LLL); 

Npts = length(envelope); tms = Dt*[0:Npts-1].';  E_t = real(envelope.*exp(-1i*E0*tms));
figure;  
subplot(2,3,[1 2 3]) 
plot(tms/T_R,E_t); xlim([tms(1),tms(end)]/T_R); xlabel('time (normalized to the round trip time)');  ylabel('Electric field');


%field autocorrelation with a 10 RT window kernel.
[F_corr,lagsF] = xcorr(E_t,E_t,Ntau); F_corr = flipud(F_corr);

display('Fourier transforming the E-field . ...'); 
NFFT = NFFTfunc(length(F_corr));win = winfunc(length(F_corr)); normconst = normalizefourier(NFFT,win,@fft);
display(['--> Non-zeropadded (decimated) field transform resolution (MHz): ' num2str(1/Dt/length(F_corr)*1E6)]);
E_w = transform(F_corr.*win,NFFT)/normconst; f_E  = 1/NFFT*[0:NFFT-1]*1/Dt; fElim = f_E > 1 & f_E < 6;

subplot(2,3,[4 5]); 
plot(f_E(fElim),abs(E_w(fElim))); xlim([3.4 4.4])
xlabel('Freq. (THz)'); ylabel(' |E(\omega)|^2 (a.u.)');

display('Fourier transforming P(t) ...'); 
P_t = abs(envelope).^2;
Y =  transform(P_t); fP = 1/Npts*[0:Npts-1]*1/Dt; 
NFFT = NFFTfunc(length(P_t));win = winfunc(length(P_t)); 
normconst = normalizefourier(NFFT,win,transform);

display(['--> Non-zeropadded P(t) transform resolution (MHz): ' num2str(1/Dt/length(P_t)*1E6)]);
P_w = transform(P_t.*win,NFFT)/normconst; 
f_P  = 1/NFFT*[0:NFFT-1]*1/Dt; 
filt = f_P >f_R/2 & f_P < 3*f_R/2;
[Ppks, Plocs,Pwidths] = findpeaks(abs(P_w(filt)),f_P(filt),'MinPeakDistance',f_R*0.75);


ax = subplot(2,3,6) ;
semilogy(f_P(filt)*1E3,abs(P_w(filt))/max(abs(P_w(filt)))); xlabel('Freq. (GHz)'); ylabel('Beatnote (a.u.)');
set(ax,'YAxisLocation','right')
xlim([(Plocs*1E3-0.5) (Plocs*1E3+0.5)]);
ylim([1E-5,1.1])






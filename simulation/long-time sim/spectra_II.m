%% clean clear and close stuff
close all; clc; 

convention = 'physics';

%determine the transform function according to the assumed convention 
if(strcmp(convention,'physics'))
    transform = @ifft;
    itransform = @fft; 
    cs =1; 
else
    transform = @fft;
    itrasform = @ifft;
    cs = -1; 
end

apodizefunc = @hanning; 
NFFTfunc = @(elem) 2^nextpow2(8*elem); 

% apodizefunc = @(elem) ones(elem,1); 
% NFFTfunc = @(elem) elem; 

%% Prep data:
%%% get the time domain envelope from the desired rtrips.
rt_start = 500; rt_end =995;

Dt = dt*iterperrecord;  iter_per_rt = round(T_R/Dt);
A_t = E_p(rt_start*iter_per_rt+1:rt_end*iter_per_rt)+E_m(rt_start*iter_per_rt+1:rt_end*iter_per_rt);
% A_t = E_p((rt_start*iter_per_rt+1:rt_end*iter_per_rt));

A_t = reshape(A_t,[length(A_t),1]); Npts = length(A_t); tms = Dt*[0:Npts-1].';  E_t = real(A_t.*exp(-cs*1i*E0*tms));




%max frequency after decimation. 

fs_dwn = 20; %in THz 
fs = 1/Dt; 
sample_freq = round(fs/fs_dwn); 

%downsample to 20 THz BW
E_t = E_t(1:sample_freq:end); tms = tms(1:sample_freq:end);
Ntau = 3*round(T_R/Dt/sample_freq);
[F_corr,lagsF] = xcorr(E_t,E_t,Ntau); F_corr = flipud(F_corr);

NFFT = NFFTfunc(length(F_corr));win = apodizefunc(length(F_corr)); normconst = normalizefourier(NFFT,win,@fft);
fprintf('Non-zeropadded field transform resolution (MHz): %.3E\n',1/Dt/sample_freq/length(F_corr)*1E6);
E_w = transform(F_corr.*win,NFFT)/normconst; f_E  = 1/NFFT*[0:NFFT-1]*1/Dt/sample_freq; fElim = f_E > 1 & f_E < 6;
dfigure; 
subplot(1,3,[1 2]); 
plot(f_E(fElim),abs(E_w(fElim))); xlim([3.4 4.4])
xlabel('Freq. (THz)'); ylabel(' |E(\omega)|^2 (a.u.)');

% subplot(2,2,2); 
% plot(lagsF*Dt*sample_freq,F_corr); 
% xlabel('Lag (ps)'); ylabel('Field autocorrelation (a.u.)');


P_t = abs(A_t).^2;
Y =  transform(P_t); fP = 1/Npts*[0:Npts-1]*1/dt; 
NFFT = NFFTfunc(length(P_t));win = apodizefunc(length(P_t)); 
normconst = normalizefourier(NFFT,win,transform);
fprintf('Non-zeropadded power transform resolution (MHz): %.3E\n',1/Dt/length(P_t)*1E6);
P_w = transform(P_t.*win,NFFT)/normconst; 
f_P  = 1/NFFT*[0:NFFT-1]*1/Dt; 
filt = f_P >f_R/2 & f_P < 3*f_R/2;
[Ppks, Plocs,Pwidths] = findpeaks(abs(P_w(filt)),f_P(filt),'MinPeakDistance',f_R*0.75);


ax = subplot(1,3,3) ;
semilogy(f_P(filt)*1E3,abs(P_w(filt))/max(abs(P_w(filt)))); xlabel('Freq. (GHz)'); ylabel('Beatnote (a.u.)');
set(ax,'YAxisLocation','right')
xlim([(Plocs*1E3-0.5) (Plocs*1E3+0.5)]);
ylim([1E-5,1.1])






%% clean clear and close stuff
close all; clc; clear;

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
NFFTfunc = @(elem) 2^nextpow2(4*elem); 

% apodizefunc = @(elem) ones(elem,1); 
% NFFTfunc = @(elem) elem; 


iterperrecord = 1;
%% Prep data:
load('CHCKPT_OPTICA(TL-sims)_(CURRENTTEST-OPTICA)_N_TRANSMISSION_LINE_4000_FP','record_U','record_V','record_v_TL','record_i_TL','record_J_TL','dat','settings'); 
E0 =dat.E0;dt =dat.dt;T_R = dat.T_R; f_R = dat.f_R;
display(['current is: ' num2str(dat.i0)])
display(['modulation amplitude is: ' num2str(settings.modA) ]);
display(['modulation freq is: ' num2str(settings.modF*dat.f_R*1E3) ]);
display(['nRF is: ' num2str(settings.nRF) ]);


%%% get the time domain envelope from the desired rtrips.
rt_start = 80; rt_end =99;

Dt = dat.dt*iterperrecord;  iter_per_rt = round(dat.T_R/Dt);
A_t = record_U(rt_start*iter_per_rt+1:rt_end*iter_per_rt)+record_V(rt_start*iter_per_rt+1:rt_end*iter_per_rt);
% A_t = E_p((rt_start*iter_per_rt+1:rt_end*iter_per_rt));
voltage = record_v_TL(rt_start*iter_per_rt+1:rt_end*iter_per_rt);

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
subplot(2,3,[1 2]);
plot(tms/T_R,abs(E_t).^2); xlim([tms(1) tms(end)]/T_R);xlabel('times (normalized to T_{Rt})'); ylabel('E(t)');

subplot(2,3,[3]);
plot(tms/T_R,voltage(1:sample_freq:end)); xlim([tms(1) tms(end)]/T_R);xlabel('times (normalized to T_{Rt})'); ylabel('V(t)');


subplot(2,3,[4 5]); 
plot(f_E(fElim),abs(E_w(fElim))); xlim([2 4])
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
display(['beatnote peak detected at: ' num2str(Plocs*1E3) 'GHz']);
display(['Full width at half prominence width is: ' num2str(Pwidths*1E6) 'MHz']);


ax = subplot(2,3,6) ;
semilogy(f_P*1E3,abs(P_w)/max(abs(P_w))); xlabel('Freq. (GHz)'); ylabel('Beatnote (a.u.)');
set(ax,'YAxisLocation','right')
xlim([1 20]);
ylim([1E-5,1.1])






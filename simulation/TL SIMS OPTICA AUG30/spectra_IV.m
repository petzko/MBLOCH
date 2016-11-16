%% clean clear and close stuff
% close all;
clc; clear;
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
load('OPTICA(GHOSTCELL-newDRY-petzko)_(12-3-3-strong)_N_3000_FP_1000','record_U','record_V','record_v_TL','record_i_TL','record_J_TL','dat','settings'); 
E0 =dat.E0;dt =dat.dt;T_R = dat.T_R; f_R = dat.f_R;
display(['central voltage is: ' num2str(settings.voltage)])
display(['modulation amplitude is: ' num2str(settings.modA) ]);
display(['modulation freq is: ' num2str(settings.modF) 'x' num2str(dat.f_R*1E3) '=' num2str(settings.modF*dat.f_R*1E3) ]);
display(['nRF is: ' num2str(settings.nRF) ]);

%%% get the time domain envelope from the desired rtrips.
rt_start = 500; rt_end = 999;

Dt = dat.dt*iterperrecord;  iter_per_rt = round(dat.T_R/Dt);
A_t = record_U(rt_start*iter_per_rt+1:rt_end*iter_per_rt)+record_V(rt_start*iter_per_rt+1:rt_end*iter_per_rt);
% A_t = E_p((rt_start*iter_per_rt+1:rt_end*iter_per_rt));
voltage = record_v_TL(rt_start*iter_per_rt+1:rt_end*iter_per_rt);
current = record_i_TL(rt_start*iter_per_rt+1:rt_end*iter_per_rt);
vs = current*dat.Z_kV_A*dat.width_mm+voltage*dat.Lp_mm;

A_t = reshape(A_t,[length(A_t),1]); Npts = length(A_t); tms = Dt*[0:Npts-1].';  E_t = real(A_t.*exp(-cs*1i*E0*tms));

%max frequency after decimation. 
fs_dwn = 20; %in THz 
fs = 1/Dt; 
sample_freq = 1 %round(fs/fs_dwn); 

%downsample to 20 THz BW
E_t = E_t(1:sample_freq:end); tms = tms(1:sample_freq:end);

NFFT = NFFTfunc(length(E_t));win = apodizefunc(length(E_t)); normconst = normalizefourier(NFFT,win,@fft);
fprintf('Non-zeropadded field transform resolution (MHz): %.3E\n',1/Dt/sample_freq/length(E_t)*1E6);
E_w = transform(E_t.*win,NFFT)/normconst; f_E  = 1/NFFT*[0:NFFT-1]*1/Dt/sample_freq; fElim = f_E > 1 & f_E < 6;

dfigure;
subplot(3,1,[1]);
plot(tms/T_R,abs(E_t).^2); xlim([tms(1) tms(end)]/T_R);xlabel('times (normalized to T_{Rt})'); ylabel('I(t)');

subplot(3,1,2); 
plot(tms/T_R,current(1:sample_freq:end)); xlim([tms(1) tms(end)]/T_R);xlabel('times (normalized to T_{Rt})'); ylabel('i_0(t)')


subplot(3,1,[3]);
plot(tms/T_R,voltage(1:sample_freq:end)); xlim([tms(1) tms(end)]/T_R);xlabel('times (normalized to T_{Rt})'); ylabel('V(t)');
dfigure;
subplot(1,2,1);
%     plot(tms/T_R,abs(E_t).^2); xlim([tms(1) tms(end)]/T_R);
%     plot(f_E(fElim(1:4:end)),abs(E_w(fElim(1:4:end))).^2); xlabel('freq (THz)'); ylabel(['|E(\omega)|^2 @' num2str(fid)]);
plot(f_E(1:20:end),abs(E_w(1:20:end))); xlim([1.5 3]); xlabel('freq (THz)'); ylabel(['|E(\omega)|^']);
     
P_t = abs(A_t).^2;
NFFT = NFFTfunc(length(P_t));win = apodizefunc(length(P_t));
normconst = normalizefourier(NFFT,win,transform);
fprintf('Non-zeropadded power transform resolution (MHz): %.3E\n',1/Dt/length(P_t)*1E6);
P_w = transform(P_t.*win,NFFT)/normconst;
V_w = transform(voltage.*win,NFFT)/normconst;
f_P  = 1/NFFT*[0:NFFT-1]*1/Dt;
filt = f_P >f_R/2 & f_P < 5*f_R/2;
[Ppks, Plocs,Pwidths] = findpeaks(abs(P_w(filt)),f_P(filt),'MinPeakDistance',f_R*0.75);
display(['beatnote peak detected at: ' num2str(Plocs*1E3) 'GHz']);
display(['Full width at half prominence width is: ' num2str(Pwidths*1E6) 'MHz']);

subplot(1,2,2);

ax = plotyy(f_P(filt)*1E3,abs(P_w(filt)),f_P(filt)*1e3,abs(V_w(filt))); xlabel('Freq. (GHz)'); ylabel(['Beatnote ' ]);
set(ax(1),'XLim',[12 15])
set(ax(2),'XLim',[12 15])

getframe
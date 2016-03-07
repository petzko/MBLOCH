%% clean clear and close stuff
close all; clc; 

convention = 'physics';

%determine the transform function according to the assumed convention 
if(strcmp(convention,'physics'))
    transform = @ifft;
    cs =1; 
else
    transform = @fft;
    cs = -1; 
end

apodizefunc = @hanning; 
NFFTfunc = @(elem) 2^nextpow2(2*elem); 

% apodizefunc = @(elem) ones(elem,1); 
% NFFTfunc = @(elem) elem; 

%% Prep data:
%%% get the time domain envelope from the desired rtrips.
rt_start = 50; rt_end =95;

Dt = dt*iterperrecord;
iter_per_rt = round(T_R/Dt);
A_t = E_p(rt_start*iter_per_rt+1:rt_end*iter_per_rt)+E_m(rt_start*iter_per_rt+1:rt_end*iter_per_rt);
% A_t = E_p((rt_start*iter_per_rt+1:rt_end*iter_per_rt));

A_t = reshape(A_t,[length(A_t),1]);
Npts = length(A_t); tms = Dt*[0:Npts-1].';
E_t = real(A_t.*exp(-cs*1i*E0*tms));

% %% Part I: calculate and plot the field spectrum.
% Nr of fft points, apodizing function of choice (hamming, box car ...) transform resolution and transform bandwidth
NFFT = NFFTfunc(Npts); win = apodizefunc(Npts); nconst = normalizefourier(NFFT,win,transform); 

% Field spectrum and the corresponding frequencies
E_w = transform(E_t.*win,NFFT)/(nconst/2);
fs = 1/Dt;
f_w  = fs/NFFT*[0:NFFT/2-1,-NFFT/2:-1].'; df = f_w(2)-f_w(1);
f_w = f_w(1:NFFT/2); E_w = E_w(1:NFFT/2); 


fprintf('Field transform bandwidth (THz): %.3E\n',fs/2);  fprintf('Non-zeropadded field transform resolution (MHz): %.3E\n',1/dt/Npts*1E6);

% power spectrum and the corresponding RF frequencies
f_R = 1/T_R;
P_t = abs(A_t).^2; 
Y_note =transform(P_t.*win,NFFT)/(nconst/2); f_note = fs/NFFT*[0:NFFT/2-1,-NFFT/2:-1].';
Y_note = Y_note(1:NFFT/2); f_note = f_note(1:NFFT/2); df = f_note(2)-f_note(1);
filt = f_note >f_R/2;
[Ppks, Plocs,Pwidths] = findpeaks(abs(Y_note(filt)),f_note(filt),'MinPeakDistance',f_R/2);
fprintf('Beatnote transform badwidth (THz): %.3E\n',fs/2); fprintf('Beatnote transform resolution (MHz): %.3E\n',1/dt/Npts*1E6);
%%% finally, do the plotting
dfigure('DName','Detected RF spectrum'); semilogy(f_note/f_R,abs(Y_note)); xlim([.5,5]);
text(Plocs(1:5)/f_R+0.02,Ppks(1:5),[num2str(Pwidths(1:5)*1E6)]);
xlabel('Freq. (per RT freq)'); ylabel('|P(w)|');  

dfigure('DName','optical field power spectrum & phase');
[Epks,Elocs,Ewidths] = findpeaks(abs(E_w/max(abs(E_w))).^2,f_w,'MinPeakHeight',0.001); 
plot(f_w,abs(E_w).^2); title('|E(\omega)|^2');  xlabel('Freq. (THz)'); ylabel('|E(w)|^2');  xlim([2,5]);
% text(Elocs,Epks*max(abs(E_w)),num2str(Ewidths*1E6));
dfigure; plot(tms/T_R,(E_t)); title('field');

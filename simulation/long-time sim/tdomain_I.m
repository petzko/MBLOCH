%% clean, clear and close stuff
close all;  clc;

convention = 'physics';

%determine the transform function according to the assumed convention 
if(strcmp(convention,'physics'))
    transform = @ifft;
    cs = 1; 
else
    transform = @fft;
    cs =-1 ; 
end

apodizefunc = @hanning; 
NFFTfunc = @(elem) 1E6; 




%% Prep data:

%%% get the time domain envelope in the desired rtrips. (max is 1000 rtrips)
rt_start =980; rt_end = 995;
Dt=  dt*iterperrecord;
iter_per_rt = round(T_R/Dt);
% A_t = E_p(rt_start*iter_per_rt+1:rt_end*iter_per_rt) + E_m(rt_start*iter_per_rt+1:rt_end*iter_per_rt);
A_t = E_p((rt_start*iter_per_rt+1:rt_end*iter_per_rt));
A_t = reshape(A_t,length(A_t),1 ); 

%total nr of sampling points,sampling time, sampling time array
N_A = length(A_t); T_tot = N_A*Dt;

% electric field at z = 0; (f_0 is our carrier frequency).
times = T_tot*linspace(-1/2,1/2,N_A)'; f0 = E0/2/pi;
E_t = real(A_t.*exp(-cs*1i*E0*times)); 

%% I Short time Fourier 
%window size of "x" number of iterations and duration in [ps]
winsize = 500; winduration = winsize*Dt;
% total number of windows
numwin = floor(N_A/winsize);

%increase the fft resolution a bit...
NFFT = NFFTfunc(winsize); win = apodizefunc(winsize); nconst = normalizefourier(NFFT,win,transform); 

t_0 = -250; t_1= 250; 
%array to store the spectra for each window.
% SPECTRA = zeros(NFFT/2,numwin);
% for w = (0:numwin-1)
%     % get the  electric field and window it with zero padding and a hanning apodization function.
%     sig = E_t(w*winsize+1:(w+1)*winsize);
%     Y = transform(sig.*win,NFFT)/(nconst/2);
%     Y = Y(1:NFFT/2);
%     SPECTRA(:,w+1) = Y;
% end
% 
% x_axis = winduration*numwin*linspace(-1/2,1/2,numwin); %in ps
% y_axis = 1/Dt*linspace(0,1/2,NFFT/2); % in THz


% figure;
% subplot(2,1,1);
% con_ax = contourf(x_axis,y_axis,abs(SPECTRA).^2);
% xlim([t_0,t_1]); xlabel(' Time (ps)'); ylim([3.4,4.4]) , ylabel('Freq. (THz)');

%% Part II
%in this part plot the time-domain intensities of the corresponding high and low
%frequency lobes. The cutoff frequency is chosen as the central freq f_0.


% Field spectrum and the corresponding frequencies
NFFT = NFFTfunc(N_A); win = apodizefunc(N_A); nconst = normalizefourier(NFFT,win,transform); 
E_w = transform(E_t)/(nconst/2); E_w = E_w(1:length(E_w)/2); 
fs = 1/Dt;  f_w = fs*linspace(0,1/2,length(E_w));


% filters
Fnorm = f0/(fs/2); % Normalized frequency
%%apply filters and correct for phase delay 
forder = 50; 
df_1 = fir1(forder,[0.00001,Fnorm]);
df_2 = fir1(forder,[Fnorm,0.99]);
E_1 = filter(df_1,1,E_t); E_2 = filter(df_2,1,E_t); 

delay = mean(grpdelay(df_1));
E_t_delayed = E_t(1:end-delay); 
tm = times(1:end-delay); 
E_1(1:delay) = []; E_2(1:delay) = []; 
E_1_env = (hilbert(E_1)); E_2_env = (hilbert(E_2));

dfigure;
subplot(2,1,2); 
plot(tm,abs(E_1_env).^2,tm,abs(E_2_env).^2); xlim([t_0,t_1]); legend('Low freq lobe','High freq lobe'); xlabel(' Time (ps)');   ylabel('Intensity (a.u.)');

NFFT= NFFTfunc(length(E_1)); win = apodizefunc(length(E_1)); nconst = normalizefourier(NFFT,win,transform); 
Y_1 = transform(E_1)/(nconst/2); Y_1 = Y_1(1:length(Y_1)/2); Y_2 = transform(E_2)/(nconst/2); Y_2 = Y_2(1:length(Y_2)/2); 
f_w = fs*linspace(0,1/2,length(Y_2)); norm = max(abs(Y_2));

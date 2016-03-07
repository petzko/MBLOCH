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
rt_start =1980; rt_end = 1995; Dt=  dt*iterperrecord; iter_per_rt = round(T_R/Dt);

% A_t = E_p(rt_start*iter_per_rt+1:rt_end*iter_per_rt) + E_m(rt_start*iter_per_rt+1:rt_end*iter_per_rt);
A_t = E_p((rt_start*iter_per_rt+1:rt_end*iter_per_rt)); A_t = reshape(A_t,length(A_t),1 ); 

%total nr of sampling points,sampling time, sampling time array
N_A = length(A_t); 

% electric field at z = 0; (f_0 is our carrier frequency).
tms = Dt*[0:N_A-1].'; f0 = E0/2/pi;
E_t = real(A_t.*exp(-cs*1i*E0*tms)); 


%% Part II
%in this part plot the time-domain intensities of the corresponding high and low
%frequency lobes. The cutoff frequency is chosen as the central freq f_0.


% Field spectrum and the corresponding frequencies
NFFT = NFFTfunc(N_A); win = apodizefunc(N_A); nconst = normalizefourier(NFFT,win,transform); 
E_w = transform(E_t,NFFT)/(nconst/2); E_w = E_w(1:length(E_w)/2);  fs = 1/Dt;  f_w = fs/NFFT*[0:NFFT/2-1].';

% filters
Fnorm = f0/(fs/2); % Normalized frequency
%%apply filters and correct for phase delay 
forder = 50; 
df_1 = fir1(forder,[0.00001,Fnorm]);
df_2 = fir1(forder,[Fnorm,0.99]);
E_1 = filter(df_1,1,E_t); E_2 = filter(df_2,1,E_t); 

delay = mean(grpdelay(df_1));
E_t_delayed = E_t(1:end-delay); 
tm = tms(1:end-delay); 
E_1(1:delay) = []; E_2(1:delay) = []; 
E_1_env = (hilbert(E_1)); E_2_env = (hilbert(E_2));

dfigure; 
subplot(2,1,1);
plot(tm/T_R,abs(E_1_env).^2,tm/T_R,abs(E_2_env).^2); legend('Low freq lobe','High freq lobe'); xlabel(' Time (in units of T_{rt})');   ylabel('Intensity (a.u.)');
xlim([5 7])


M = 50; 

E_a = hilbert(E_t);
E_r = real(E_a); E_i = imag(E_a); 
% 
% 
f1 = 1/(2*pi*dt)*atan((E_r(1:end-1).*E_i(2:end)-E_r(2:end).*E_i(1:end-1))./(E_r(1:end-1).*E_r(2:end)+E_i(1:end-1).*E_i(2:end)));
% f1_avg = mvavg(f1,M); 
f2 = 1/(4*pi*dt)*atan((E_r(1:end-2).*E_i(3:end)-E_r(3:end).*E_i(1:end-2))./(E_r(1:end-2).*E_r(3:end)+E_i(1:end-2).*E_i(3:end)));
% f2_avg = mvavg(f2,M);

subplot(2,1,2);
% plot(tms(2:end-1)/T_R,f2_avg,'LineWidth',0.7);
hold on;  
plot(tms(1:end-1)/T_R,f1,tms(2:end-1)/T_R,f2); 
hold off; xlabel('Time (in units of T_{rt})'); ylabel('Frequency (THz)'); xlim([5 7]);ylim([3 5]);


% 
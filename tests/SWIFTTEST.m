 
clear; close;
%periodic pulse of gaussians with repetition rate f_r (period T = 1/f_r); 

%repetition rate/Rt Freq!
fr = 0.45; % 2 times per sec;
Tr = 1/fr;

FWHM = 0.3; %sec
sig = FWHM/(2*sqrt(2*log(2))); % std. dev

Tend = 100*Tr ; 
N = 2^18; 

centers = 0:Tr:Tend;
gauss = @(t,mu) exp(-(t-mu).^2/(2*sig^2))/(sig*sqrt(2*pi));
times = linspace(0,Tend,N)';
dt = times(2)-times(1);

signal = zeros(N,1);

for i = 1:length(centers)
    signal = signal + gauss(times,centers(i));
end

% signal = shiftRealSig(signal,dt,1,30);
window_duration = 10*Tr; window_vector_len = window_duration/dt;

num_win = floor(length(signal)/window_vector_len); 
options = {'no','truncate'};

sig_windowed = signal(1:window_vector_len-1);
sig_avg = sig_windowed;
dv = fr;

[ SP_GLOB,freq_glob,Y01,Y02 ] = SpectrumProduct( sig_windowed,dt,1,dv,1,options);
[SWIFT,S_I,S_Q,lags] = SWIFTinterferogramFast(sig_windowed,dt,1,dv,options);
[freqswift_glob,YSWIFT_GLOB ] = FourierTransf((SWIFT),dt,0,1,1,options);

coherence = abs(YSWIFT_GLOB)./SP_GLOB;
apodize = hanning(window_vector_len-1);

for n = 2:(num_win)
    sig_windowed = signal((n-1)*window_vector_len+1:n*window_vector_len-1);
    sig_avg = sig_avg+sig_windowed;
    [SP,muddy1,muddy2,muddy3] = SpectrumProduct(sig_windowed,dt,1,dv,1,options);
    [SWIFT,dummy1,dummy2,dummy3] = SWIFTinterferogramFast(sig_windowed,dt,1,dv,options);
    [dummy,Yswift ] = FourierTransf((SWIFT),dt,0,1,1,options);
    SP_GLOB = SP_GLOB+SP;
    coherence  = coherence+abs(Yswift)./SP;
    YSWIFT_GLOB = YSWIFT_GLOB + Yswift; 
end
    sig_avg = sig_avg/num_win; YSWIFT_GLOB = YSWIFT_GLOB/num_win; SP_GLOB = SP_GLOB/num_win;
% calculate the coherece! 

coherence = 0*YSWIFT_GLOB; partition_size = floor(length(freqswift_glob)*(3*fr)/freqswift_glob(end)); 
num_partitions = floor(length(freqswift_glob)/partition_size);
for i =1:num_partitions
    interval = (i-1)*partition_size+1:partition_size*i;
    coherence(interval) = abs(YSWIFT_GLOB(interval))./SP_GLOB(interval);
end

% coherence = coherence.*(coherence <=1);

%%% swift analysis
subplot(2,1,1); 
semilogy(freq_glob,(SP_GLOB),freqswift_glob,abs(YSWIFT_GLOB));

subplot(2,1,2); 
P_t = abs(signal).^2; %sum and difference frequencies! 
[f_w, P_w] = FourierTransf(P_t,dt,0,1,1,options);
plot(f_w,abs(P_w));

figure; semilogy(freq_glob,coherence);





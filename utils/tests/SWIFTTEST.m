

%periodic pulse of gaussians with repetition rate f_r (period T = 1/f_r); 

%repetition rate
fr = 0.5; % 2 times per sec;
Tr = 1/fr; % rt freq! 
f0 = 0; %hz

FWHM = 0.3; %sec
sig = FWHM/(2*sqrt(2*log(2))); % std. dev

Tend = 500 ; 
N = 2^18; 

centers = 0:Tr:Tend;
gauss = @(t,mu) exp(-(t-mu).^2/(2*sig^2))/(sig*sqrt(2*pi));
times = linspace(0,Tend,N)';
dt = times(2)-times(1);

signal = zeros(N,1);

for i = 1:length(centers)
    signal = signal + gauss(times,centers(i));
end

signal = shiftRealSig(signal,dt,1,100);

subplot(4,1,1);
plot(times,(signal));

title('signal');
options = {'hamming', 'truncate'};

[ SP,freq,Y01,Y02 ] = SpectrumProduct( signal,dt,1,fr,1,options);
[SWIFT,S_I,S_Q,lags] = SWIFTinterferogramFast(signal,dt,1,0.98*fr,options);
[freqswift,Yswift ] = FourierTransf((SWIFT),dt,f0,1,1,options);



subplot(4,1,2); 
plot(freq,abs(SP));
title('spectrum product');

subplot(4,1,3); 
plot(freqswift,abs(Yswift));
title('swift');
subplot(4,1,4); 
semilogy(freq,abs(SP),freqswift,abs(Yswift));
title('SP & SWIFT');






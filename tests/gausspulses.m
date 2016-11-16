 
clear; close;
%periodic pulse of gaussians with repetition rate f_r (period T = 1/f_r); 

%repetition rate/Rt Freq!
fr = 0.5; % 2 times per sec;
Tr = 15/fr;

FWHM = 5; %sec
sig = FWHM/(2*sqrt(2*log(2))); % std. dev

Tend = 100*Tr ; 
N = 2^18; 

freq = .5;

centers = 0:Tr:Tend;
gauss = @(t,mu) ...
    exp(-(t-mu).^2/(2*sig^2))/(sig*sqrt(2*pi)).*exp(1i*freq*2*pi*t);
times = linspace(0,Tend,N)';
dt = times(2)-times(1);

signal = zeros(N,1);

for i = 1:length(centers)
    signal = signal + gauss(times,centers(i));
end
subplot(2,1,1)
plot(times,real(signal),times,abs(signal))
Y1 =  fft(signal);

f_rep = 1/Tr;
fc = freq;
m =.5

fm_signal = @(t) cos(2*pi*fc*t+10*cos(2*pi*f_rep*t));
signal = fm_signal(times);
subplot(2,1,2);
plot(times,real(signal))
Y2 = fft(signal);

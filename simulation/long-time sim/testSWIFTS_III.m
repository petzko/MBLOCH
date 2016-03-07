
clear; close;
%periodic pulse of gaussians with repetition rate f_r (period T = 1/f_r); 

%repetition rate/Rt Freq!
fr = 0.45; % 2 times per sec;
Tr = 1/fr;

FWHM = 0.3; %sec
sig = FWHM/(2*sqrt(2*log(2))); % std. dev

Tend = 100*Tr ; 
N = 2^20; 

centers = 0:Tr:Tend;
gauss = @(t,mu) exp(-(t-mu).^2/(2*sig^2))/(sig*sqrt(2*pi));
times = linspace(0,Tend,N)';
dt = times(2)-times(1);
signal = zeros(N,1);
for i = 1:length(centers)
    signal = signal + gauss(times,centers(i));
end

E0 = 30*2*pi;


params.N = length(signal); 
params.dt = dt; 
params.iterperrec = 1; 
params.dv = fr;
params.E0 = E0;
params.complex = 'no'
params.convention = 'physics'; 
iter_per_rt = round(Tr/dt);
Ntau = round(20*iter_per_rt); params.Ntau = Ntau; 

data = calcswift( signal, params );
% subplot(3,1,1);  plot(data.lags0,data.S0);
% subplot(3,1,2);  plot(data.lagsI, data.SI);
% subplot(3,1,3);  plot(data.lagsQ, data.SQ);
plot(data.f0,abs(data.P0),data.f0,abs(data.P1));

%% plot swifts 
figure; 
semilogy(data.f0,sqrt(abs(data.P0).*abs(data.P1)),data.f0,abs(data.SWIFT));

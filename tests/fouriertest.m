dt = 1e-2;
T = 1.5;
t = 0:dt:T;
f0 = 0;
fs = 1/dt;
f_1 = 3;
df = 5;
f_2 = 9;

xcomplex  = 3*exp(1i*(2*pi*f_1*t + 0.1)) + 4*exp(1i*(2*pi*f_2*t + 0.25));
xreal = 3*sin(2*pi*f_1*t) + 4*sin(2*pi*f_2*t);

xcos = dct(xreal);
xsin = dst(xreal);

xnew = shift(xcos,[0 3],0);



xreal_new = idct(xnew); 


plot(xreal); 
hold on; 
plot(xreal_new);

[freqR,YR] = FourierTransf(xreal,dt,f0,1,1);
[freqR_new,YR_new] = FourierTransf(xreal_new,dt,f0,1,1);
figure; 
hold on;
plot(freqR,abs(YR));
plot(freqR_new,abs(YR_new));


% 
% [freqC,YC] = FourierTransf(xcomplex,dt,f0,1,1);

% 
% subplot(2,1,1);
% plot(freqR,abs(YR));
% title('real');
% subplot(2,1,2);
% plot(freqC,abs(YC));
% title('complex');
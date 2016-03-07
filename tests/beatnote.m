clear;clc;close all;
f_11 = 205; A_11 = 1;
f_12 = 204; A_12 = 2; 

f_21 = 104; A_21 = 1;
f_22 = 102; A_22 = 5;

s_1 = @(t) A_11 * exp(1i*2*pi*f_11*t) + A_12 *exp(1i*2*pi*f_12*t) ;
s_2 = @(t) A_21 * exp(1i*2*pi*f_21*t) + A_22 * exp(1i*2*pi*f_22*t);

s_beat = @ (t ) s_1(t)+s_2(t);

t_end = 6;
t = 0; 
dt = 1E-3;
S_1 = 0; S_2=0 ;S_3 = 0;S_cos = 0;
ctr = 0;

while t< t_end
    ctr = ctr + 1;
    S_1(ctr) = s_1(t);
    S_2(ctr) = s_2(t);
    S_beat(ctr) = s_beat(t);
%     S_cos(ctr) = 2*cos(2*pi*(f_11-f_21)*t/2);
    t = t+dt;
end

%%
    times = linspace(0,t,ctr);
    subplot(3,1,1);
    plot(times,real(S_1));
    subplot(3,1,2);
    plot(times,real(S_2));
    subplot(3,1,3);
    plot(times,real(S_beat),times,S_cos,'ro');
    getframe;
%     
    N = length(S_beat); 
    apodization = hamming(N)'; 
    NFFT =  2^nextpow2(N);
    Y_1 = fft((S_beat),NFFT); Y = Y_1(1:length(Y_1)/2)/length(Y_1); 
    
    f = 1/dt*linspace(0,1/2,length(Y));
    figure; 
    plot(f,abs(Y));
    

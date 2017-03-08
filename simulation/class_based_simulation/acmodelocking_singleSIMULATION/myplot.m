clc;clear;close all;
% load variables from mat file
name1 = matfile('name1.mat');
name2 = matfile('name2.mat');
name3 = matfile('name3.mat');
name4 = matfile('name4.mat');
name5 = matfile('name5.mat');
name6 = matfile('name6.mat');



dt = name1.dt; % unit: ps
T_R = name1.T_R; % round trip time, unit: ps
t1 = round(150*T_R/dt);% set plot time range (how many round trips)
t2 = round(155*T_R/dt);


% import electric feld vector
record_U1 = name1.record_U;
record_U2 = name2.record_U;
record_U3 = name3.record_U;
record_U4 = name4.record_U;
record_U5 = name5.record_U;
record_U6 = name6.record_U;

% import voltage vector
record_v_TL1 = name1.record_v_TL;
record_v_TL2 = name2.record_v_TL;
record_v_TL3 = name3.record_v_TL;
record_v_TL4 = name4.record_v_TL;
record_v_TL5 = name5.record_v_TL;
record_v_TL6 = name6.record_v_TL;


%******** plot pulse --> electric feld and bias in time domain ****************
figure
x = dt*linspace(t1,t2,t2-t1+1)/T_R;
E1 = abs(record_U1(t1:t2));
subplot(6,1,1)
plot(x,E1)
ylim([0,0.7])
title('Spectra in time domain (5 round trips with t_{rt} = 72.0498 ps, f_{RF} = f_{rt})')
ylabel('|E_z|')
legend('modA = 0')

E2 = abs(record_U2(t1:t2));
subplot(6,1,2)
plot(x,E2)
ylim([0,0.7])
ylabel('|E_z|')
legend('modA = 0.02')

E3 = abs(record_U3(t1:t2));
subplot(6,1,3)
plot(x,E3)
ylim([0,0.7])
ylabel('|E_z|')
legend('modA = 0.04')

E4 = abs(record_U4(t1:t2));
subplot(6,1,4)
plot(x,E4)
ylim([0,0.7])
xlabel('round trips')
ylabel('|E_z|')
legend('modA = 0.06')

E5 = abs(record_U5(t1:t2));
subplot(6,1,5)
plot(x,E5)
ylim([0,0.7])
xlabel('round trips')
ylabel('|E_z|')
legend('modA = 0.08')

E6 = abs(record_U6(t1:t2));
subplot(6,1,6)
plot(x,E6)
ylim([0,0.7])
xlabel('round trips')
ylabel('|E_z|')
legend('modA = 0.10')


%******** plot spectra --> in frequency domain ****************
[ F1,f1 ] = mydft2(record_U1(600000:end),dt);
[ F2,~ ] = mydft2(record_U2(600000:end),dt);
[ F3,~ ] = mydft2(record_U3(600000:end),dt);
[ F4,~ ] = mydft2(record_U4(600000:end),dt);
[ F5,~ ] = mydft2(record_U5(600000:end),dt);
[ F6,~ ] = mydft2(record_U6(600000:end),dt);


f1 = f1/1E12;  % Frequency unit: GHz

figure
subplot(6,1,1)
plot(f1,normc(abs(F1)))
xlim([-1,1])
title('Spectra in frequency domain (f_{RF} = f_{rt} ) = 13.9 Ghz')
ylabel('normalized intensity')
legend('modA = 0')


subplot(6,1,2)
plot(f1,normc(abs(F2)))
xlim([-1,1])
xlabel('frequency /THz')
ylabel('normalized intensity')
legend('modA = 0.02')


subplot(6,1,3)
plot(f1,normc(abs(F3)))
xlim([-1,1])
xlabel('frequency /THz')
ylabel('normalized intensity')
legend('modA = 0.04')


subplot(6,1,4)
plot(f1,normc(abs(F4)))
xlim([-1,1])
xlabel('frequency /THz')
ylabel('normalized intensity')
legend('modA = 0.06')

subplot(6,1,5)
plot(f1,normc(abs(F5)))
xlim([-1,1])
xlabel('frequency /THz')
ylabel('normalized intensity')
legend('modA = 0.08')

subplot(6,1,6)
plot(f1,normc(abs(F6)))
xlim([-1,1])
xlabel('frequency /THz')
ylabel('normalized intensity')
legend('modA = 0.10')


%******** plot THz field beatnot and voltage over frequncy ****************
P1 = abs(record_U1(400000:end)).^2;
P2 = abs(record_U2(400000:end)).^2;
P3 = abs(record_U3(400000:end)).^2;
P4 = abs(record_U4(400000:end)).^2;
P5 = abs(record_U5(400000:end)).^2;
P6 = abs(record_U6(400000:end)).^2;

[ B1,f2 ] = mydft2(P1,dt); B1 = B1(1:length(B1)/2);
[ B2,~ ] = mydft2(P2,dt); B2 = B2(1:length(B2)/2);
[ B3,~ ] = mydft2(P3,dt); B3 = B3(1:length(B3)/2);
[ B4,~ ] = mydft2(P4,dt); B4 = B4(1:length(B4)/2);
[ B5,~ ] = mydft2(P5,dt); B5 = B5(1:length(B5)/2);
[ B6,~ ] = mydft2(P6,dt); B6 = B6(1:length(B6)/2);

% [ V1,~ ] = mydft2(record_v_TL1(600000:end)-bias*ones(length(record_v_TL1(600000:end))),dt); V1 = V1(1:length(V1)/2);
% [ V2,~ ] = mydft2(record_v_TL2(600000:end)-bias*ones(length(record_v_TL2(600000:end))),dt); V2 = V2(1:length(V2)/2);
% [ V3,~ ] = mydft2(record_v_TL3(600000:end)-bias*ones(length(record_v_TL3(600000:end))),dt); V3 = V3(1:length(V3)/2);
% [ V4,~ ] = mydft2(record_v_TL4(600000:end)-bias*ones(length(record_v_TL4(600000:end))),dt); V4 = V4(1:length(V4)/2);
tend = length(record_U1);
tx = dt*linspace(400000,tend,tend-400000+1);
ty = sin(2*pi/T_R*tx);
[ V1,~ ] = mydft2(ty,dt);   V1 = V1(1:length(V1)/2);

f2 = f2'/1e9; f2=f2(1:length(f2)/2);


figure
subplot(6,1,1)
% semilogy(f2,normc(abs(B1)))
yyaxis left
plot(f2,normc(abs(B1)))
title('THz field beatnote at mod freq 13.9 Ghz')
xlim([10,20])
yyaxis right
plot(f2,abs(V1))
% ylabel('Beatnote signal','RF signal')
legend('modA = 0')


subplot(6,1,2)
% semilogy(f2,normc(abs(B2)))
yyaxis left
plot(f2,normc(abs(B2)))
xlim([10,20])
% ylabel('Beatnote signal','RF signal')
yyaxis right
plot(f2,abs(V1))
legend('modA = 0.02')


subplot(6,1,3)
% semilogy(f2,normc(abs(B3)))
yyaxis left
plot(f2,normc(abs(B3)))
xlim([10,20])
% ylabel('Beatnote signal','RF signal')
yyaxis right
plot(f2,abs(V1))
legend('modA = 0.04')


subplot(6,1,4)
% semilogy(f2,normc(abs(B4)))
yyaxis left
plot(f2,normc(abs(B4)))
xlim([10,20])
% ylabel('Beatnote signal','RF signal')
yyaxis right
plot(f2,abs(V1))
xlabel('frequency /GHz')
legend('modA = 0.06')


subplot(6,1,5)
% semilogy(f2,normc(abs(B4)))
yyaxis left
plot(f2,normc(abs(B5)))
xlim([10,20])
% ylabel('Beatnote signal','RF signal')
yyaxis right
plot(f2,abs(V1))
xlabel('frequency /GHz')
legend('modA = 0.08')


subplot(6,1,6)
% semilogy(f2,normc(abs(B4)))
yyaxis left
plot(f2,normc(abs(B6)))
xlim([10,20])
% ylabel('Beatnote signal','RF signal')
yyaxis right
plot(f2,abs(V1))
xlabel('frequency /GHz')
legend('modA = 0.10')




%%******** plot pulse --> electric feld and bias in time domain ****************
% figure
% x = dt*linspace(t1,t2,t2-t1+1)/T_R;
% y1 = abs(record_U1(t1:t2));
% u1 = abs(record_v_TL1(t1:t2));
% subplot(4,1,1)
% plotyy(x,y1,u1)
% ylim([0,0.7])
% title('Spectra in time domain (5 round trips with t_{RT} = 72.0498 ps)')
% ylabel('|E_z|')
% legend('modA = 0')
% 
% y2 = abs(record_U2(t1:t2));
% u2 = abs(record_v_TL2(t1:t2));
% subplot(4,1,2)
% plotyy(x,y2,u2)
% ylim([0,0.7])
% ylabel('|E_z|')
% legend('modA = 0.02')
% 
% y3 = abs(record_U3(t1:t2));
% u3 = abs(record_v_TL3(t1:t2));
% subplot(4,1,3)
% plotyy(x,y3,u3)
% ylim([0,0.7])
% ylabel('|E_z|')
% legend('modA = 0.04')
% 
% y4 = abs(record_U4(t1:t2));
% u4 = abs(record_v_TL4(t1:t2));
% subplot(4,1,4)
% plotyy(x,y4,u4)
% ylim([0,0.7])
% xlabel('round trips')
% ylabel('|E_z|')
% legend('modA = 0.06')

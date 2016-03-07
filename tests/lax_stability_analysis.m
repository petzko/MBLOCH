clear;

c = 1;
l = 0;
k = 0;

dt  = 0.01;
N = 100;
beta = linspace(0,2*pi,N);
G = @(b,r) (1+k*dt + k^2*dt*dt/2)-1i*r*(1+l*dt)*sin(b)+r^2*(cos(b)-1)

r = 0.3;
dx_03 = c*dt/r;
G_03  = G(beta,r);
phase_03 = angle(G_03);


r = 0.5;
dx_05 = c*dt/r;
G_05  = G(beta,r);
phase_05 = angle(G_05);

r = 0.7;
dx_07 = c*dt/r;
G_07  = G(beta,r);
phase_07 = angle(G_07);

r = 0.9;
dx_09 = c*dt/r;
G_09  = G(beta,r);
phase_09 = angle(G_09);

r = 1.0;
dx_10 = c*dt/r;
G_10  = G(beta,r);
phase_10 = angle(G_10);

r = 1.1;
dx_11 = c*dt/r;
G_11  = G(beta,r);
phase_11 = angle(G_11);

r = 1.3;
dx_13 = c*dt/r;
G_13  = G(beta,r);
phase_13 = angle(G_13);
beta = beta(1:3:end)
figure;
plot(beta,abs(G_03(1:3:end)),':rs',beta,abs(G_05(1:3:end)),':g+',beta,abs(G_07(1:3:end)),':b^',beta,abs(G_09(1:3:end)),':y*',beta,abs(G_10(1:3:end)),':cx',beta,abs(G_11(1:3:end)),':mo',beta,abs(G_13(1:3:end)),'k')
title(['Amplification Factor Lax scheme: k = ' num2str(k) ' , c = ' num2str(c) ' , dt = ' num2str(dt) ]);
legend('r = 0.3','r = 0.5', 'r = 0.7','r = 0.9', 'r = 1.0' , 'r = 1.1', 'r = 1.3');
axis([0 2*pi 0 3])
% 
% 
% figure 
% plot(beta,phase_03/pi,'r',beta,phase_05/pi,'g',beta,phase_07/pi,'b',beta,phase_09/pi,'y',beta,phase_10/pi,'c',beta,phase_11/pi,'m',beta,phase_13/pi,'k')
% title(['Phase error Lax Scheme: l = ' num2str(l) ' , c = ' num2str(c) ' , dt = ' num2str(dt) ]);
% legend('r = 0.3','r = 0.5', 'r = 0.7','r = 0.9', 'r = 1.0' , 'r = 1.1', 'r = 1.3');
% 
% 
% figure
% plot(phase_03(2:N-1)/pi,abs(G_03(2:N-1)),'r',phase_05(2:N-1)/pi,abs(G_05(2:N-1)),'g',phase_07(2:N-1)/pi,abs(G_07(2:N-1)),'b',phase_09(2:N-1)/pi,abs(G_09(2:N-1)),'y',phase_10(2:N-1)/pi,abs(G_10(2:N-1)),'c',phase_11(2:N-1)/pi,abs(G_11(2:N-1)),'m',phase_13(2:N-1)/pi,abs(G_13(2:N-1)),'k')
% title(['Phase error Lax Scheme: l = ' num2str(l) ' , c = ' num2str(c) ' , dt = ' num2str(dt) ]);
% legend('r = 0.3','r = 0.5', 'r = 0.7','r = 0.9', 'r = 1.0' , 'r = 1.1', 'r = 1.3');

% figure
% subplot(2,1,1);
% plot(beta,phase_13/pi);
% subplot(2,1,2);
% plot(beta,abs(G_13))
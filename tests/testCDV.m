% clear;clc; close all;
% 
% c = +1;
% a = 0; b = 5;
% beta = 1i;
% N = 1024; 
% % x = (a+b)/2-(b-a)/2*cos(pi*(0:N-1)/(N-1))';
% x = linspace(a,b,N)';
% dx = x(2)-x(1); dt = 0.1*dx/abs(c);
% E_0 = 1;
% initial_conditions = zeros(N,1);
% 
% tau =0.01; %sec
% aE = @(z,t) E_0.*sech((t-(z-b/2)/c)/tau).*exp(1i*-1.5*z);
% aE_t =  @(z,t) -E_0/tau.*sech((t-(z-b/2)/c)/tau).*tanh((t-(z-b/2)/c)/tau).*exp(1i*-1*z);
% 
% 
% Ux = aE(x,0); Vx = aE_t(x,0);
% Vx(1) = 0;  Vx(end) = 0 ; Ux(1) = 0  ; Ux(end) = 0; 
% 
% 
% dk = 2*pi*dx; k = dk*(0:1:N-1);
% diagK = diag(k); 
% 
% I = eye(N,N); 
% A = [2*1i/c/beta*I (2*diagK/beta); I zeros(N,N)];
% expA = expm(A*dt); 
% 
% ctr = 1;
% t = dt;
% 
% tend = 2*(2*(b-a)/abs(c));
% options = '';
% Eoft = 0 ;
% 
% while t<tend
%     
%     Vk = fft(Vx); Uk = fft(Ux); 
%     Yk = [Vk;Uk];
%     
%     Yk = expA*Yk ; 
%     Vx = ifft(Yk(1:N)); Ux = ifft(Yk(N+1:end));
% 
%     if(mod(ctr,1)==0)
%         plot(x,abs(Ux));
%         title(['t = ' num2str(t)]);
%         getframe;
%     end
%     
%     ctr =ctr+1;
%     t=t+dt;
% end

clear;clc; close all;

c = +1; E_0 = 0.25
b = 2; beta = 2;
Nt = 1024; 
T = 2 ; 
dt= T/Nt;
t = (-Nt/2:1:Nt/2-1)'*dt;

dw = 2*pi/T; w = [0:Nt/2-1 0 -Nt/2+1:-1]'*dw;

Z = 300;
dz = 0.03; numsteps = round(Z/dz); 

tau =0.01; %sec
Ut = 0.25.*exp(-(t/5).^2); Uf = fft(Ut);

uplot = abs(Uf); zplot = 0;

L = +1/c*w; 
N = 1i*beta/2*w.^2; 


% expN = expm(N*dt); 

ctr = 1;
t = dt;

options = '';
Eoft = 0 ;

saveinterval = 200; 
for istep = 1:numsteps

    Uf = exp((L+N)*dz).*Uf;
    
    if(mod(istep,saveinterval) == 0)
        U = ifft(Uf)';
        uplot = [uplot;abs(U)];
        zplot = [zplot;istep*dz];
        display(['step = ' num2str(istep)]);
    end


end
waterfall(t, zplot, uplot); 






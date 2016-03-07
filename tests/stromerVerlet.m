clear;clc; close all;
%Velocity, domain boundaries, dispersion coeff, grid size, grid, grid spacing and time step 
c = +2; a = -10; b = 10;
beta = +3; N = 250;  x = linspace(a,b,N)'; 
dx = x(2)-x(1); dt = 0.01*dx/abs(c);

% %differentiation matrix.  <- upwind
Dx = zeros(N,N); 
% Dx(1,1) = 1/dx; Dx(1,N-1) = -1/dx;
for i = 2:N-1 Dx(i,i) = 1/dx; Dx(i,i-1) = -1/dx; end; 

% %differentiation matrix.  <- central
% Dx = zeros(N,N); 
% Dx(1,2) = 1/2/dx; Dx(1,N-1) = -1/2/dx;
% for i = 2:N-1 Dx(i,i-1) = -1/2/dx; Dx(i,i+1) = 1/2/dx; end; 
% Dx(N,2) = 1/2/dx; Dx(N,N-1) = -1/2/dx;

% % Different initial conditions 
z0 = (b+a)/20;

% % 1) A sech pulse. 
tau =1; %FWHM duration of the pulse 
k = 5; %arbitrary phase shift... 
aE = @(z,t) sech((t-(z-z0)/c)/tau).*exp(1i*k*z); 
aE_t =  @(z,t) -1/tau.*sech((t-(z-z0)/c)/tau).*tanh((t-(z-z0)/c)/tau).*exp(1i*k*z);

% % 2) A gaussian 
% tau =1; %sec
% k = 5; %arbitrary phase shift... 
% aE = @(z,t)exp(-((t-(z-z0)/c)/tau).^2).*exp(1i*k*z);
% aE_t =  @(z,t) 2*(t-(z-z0)/c)/tau.*exp(-((t-(z-z0)/c)/tau).^2).*exp(1i*k*z);

% % 3) A travelling sine wave! 
% m = 5; %order of the oscillations
% kz = 2*pi/(b-a)*m % wave number! 
% w = kz*c;
% aE = @(z,t)sin(w*t-kz*z);
% aE_t =  @(z,t) w*cos(w*t-kz*z);

U = aE(x,0);  V = aE_t(x,0); 
U(end) = U(end-1); V(end) = V(end-1); 

%mass;  damping; "elasticity" matrix!
M = 1i*beta/2; C = (1/c); D = Dx;
ctr = 1; t = dt;

tend = 2*(2*(b-a)/abs(c));
Unew = U; Vnew = V; 
Uanal = U; Vanal = V; 

while t<tend

    if(mod(ctr-1,100)==0)
        plot(x,real(U),x,real(V),x,real(aE(x,t)),x,real(aE_t(x,t)));
%         plot(x,real(Uanal),x,real(Vanal)); %<- plot the analytical expressions!  
%         legend('U^n','V^n','U(t)','V(t)');
        title(['t = ' num2str(t)]);
        xlim([a-1 b+1])
        ylim([-10 10])
        getframe;
    end
    
    Uanal(2:end) = aE(x(2:end),t); Uanal(1) = Uanal(end);
    Vanal(2:end) = aE_t(x(2:end),t); Vanal(1) = Vanal(end);
    %acceleration, U and velocity !
    A = -(C*V+D*U)/M;
    
    
    Unew = U + V*dt+0.5*(dt^2)*A;
    Vnew = (V+dt/2*A-dt/2*D*Unew/M)/(1+dt/2*C/M);
    U = Unew; V = Vnew;  
    U(end) = U(end-1); V(end) = V(end-1);
    ctr =ctr+1;     t=t+dt;
end





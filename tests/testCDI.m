clear;clc; close all;

c = +1;
a = -1; b = 1;
beta = 10;
N = 500; 
% x = (a+b)/2-(b-a)/2*cos(pi*(0:N-1)/(N-1))';
x = linspace(a,b,N)';
dx = x(2)-x(1); dt = 0.5*dx/abs(c);
E_0 = 1;
initial_conditions = zeros(N,1);
% Dx = (2/(b - a))*SpectralMSolver.cheb(N-1);
Dx = zeros(N,N);
Dx(1,1) = 1/dx; 
for i = 2:N Dx(i,i) = 1/dx; Dx(i,i-1) = -1/dx; end; 

tau =0.1; %sec

% aE = @(z,t) E_0.*sech((t-(z-b/2)/c)/tau);
% aE_t =  @(z,t) -E_0/tau.*sech((t-(z-b/2)/c)/tau).*tanh((t-(z-b/2)/c)/tau);

aE = @(z,t) E_0.*exp(-((t-(z-b/2)/c)/tau).^2)
aE_t =  @(z,t) -1/5.*2*(t-(z-b/2)/c)/tau.* E_0.*exp(-((t-(z-b/2)/c)/tau).^2);


EnumRK = aE(x,0); 
UnumRK = aE_t(x,0);
EnumRK(1) = 0;  EnumRK(end) = 0 ; UnumRK(1) = 0  ; UnumRK(end) = 0; 

EnumCN = EnumRK; EnumExp = EnumRK;
UnumCN = UnumRK; UnumExp = UnumRK;

I = eye(N,N); 
A = [2*1i/c/beta*I ((2*1i/beta)*Dx); I zeros(N,N)];

% boundary conditions Y(1) = Y(N) and Y(N+1) = Y(2N)
P = eye(size(A)); 
if( c > 0) % what numerical boundary conditions shall I implement !??? 
    P(1,:) =  P(N,:);  % for first order PBC
    P(N+1,:) = P(end,:);
else
    P(N,:) = P(1,:);
    P(2*N,:) = P(1,:); 
end

%projection matrix! <- implements boundary conditions! 
PA = P*A;

% Cranck Nicholson matrix 
M1 = eye(size(PA))-dt/2*PA; M2 = eye(size(PA))+dt/2*PA; 
M = inv(M1)*M2;

%Matrix Exponential
Aexp = expm(PA*dt); 

ctr = 1;
t = dt;

tend = 2*(2*(b-a)/abs(c));
options = '';
Eoft = 0 ;

while t<tend

    YRK = [UnumRK;EnumRK];
    YCN = [UnumCN;EnumCN];
    Yexp = [UnumExp;EnumExp];
  
% 4th order Runge-Kutta
    k1 = propGVD(t,YRK,options,PA);
%     k1(1,:) = k1(N,:); k1(N+1,:) = k1(2*N,:);
    k2 = propGVD(t+dt/2,YRK+dt/2*k1,options,PA);
%     k2(1,:) = k2(N,:); k2(N+1,:) = k2(2*N,:);
    k3 = propGVD(t+dt/2,YRK+dt/2*k2,options,PA);
%     k3(1,:) = k3(N,:); k3(N+1,:) = k3(2*N,:);
    k4 = propGVD(t+dt,YRK+dt*k3,options,PA);
%     k4(1,:) = k4(N,:); k4(N+1,:) = k1(2*N,:);
    
    YRK = YRK+dt/6*(k1+2*k2+2*k3+k4);
    UnumRK = YRK(1:N) ; EnumRK = YRK((N+1):end);
%     EnumRK(1) = 0;  UnumRK(1) = 0; UnumRK(end) = 0; 
%     EnumRK(end) = 0*EnumRK(end-1);
%     
% % Cranck-Nicholson
%     YCN = M*YCN;
%     UnumCN = YCN(1:N) ; EnumCN = YCN((N+1):end);
% % Matrix exponential
%     Yexp = Aexp*Yexp;
%     UnumExp = Yexp(1:N) ; EnumExp = Yexp((N+1):end);
  
    if(mod(ctr,100)==0)
        plot(x,real(EnumRK),x,real(EnumCN),x,real(EnumExp));
        legend('RK','CN','EXP');
        title(['t = ' num2str(t)]);
%         axis([a b -10 10]);
        getframe;
    end
    
    ctr =ctr+1;
    t=t+dt;
end





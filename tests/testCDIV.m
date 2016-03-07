clear;clc; close all;

c = 1; n = 1; 
a = 0; b = 5;
beta = 4; N = 250; 
E_0 = 1; initial_conditions = zeros(N,1);
x = linspace(a,b,N)';
dx = x(2)-x(1);
dt = dx/(abs(c));


tau =0.01; %sec

aE = @(z,t) E_0.*sech((t-(z-(b+a)/2)/c)/tau);
aE_t =  @(z,t) -E_0/tau.*sech((t-(z-(b+a)/2)/c)/tau).*tanh((t-(z-(b+a)/2)/c)/tau);

sig2 = ((b-a)/10)^2;
k = 20;
aE = @(z,t) E_0.*exp(-(t-z+(b+a)/2).^2/2/sig2).*exp(1i*k*z);
aE_t =  @(z,t)-E_0*(t-z+(b+a)/2).*exp(-(t-z+(b+a)/2).^2/2/sig2).*exp(1i*k*z) ;


U = aE(x,0); V = aE_t(x,0);
solver = RNFDSolver(N,dx,c/abs(c),abs(c), U);
Uorig = U; 

ctr = 1;
t = dt;

tend = 100*(2*(b-a)/abs(c));
U_t = 0 ;

dummy = 0*U; 

fidx = 2:N; 
bidx = 1:N-1;
Ut = U; 

I = eye(N); 
alpha = 2*1i/beta/c; 

A = [alpha*I zeros(N); I zeros(N)];
options = '';

while t<tend
    
    
    
    U = solver.make_step(dummy,dummy,dummy,dt);
    U = solver.set_bdry( U(end),'no');
  
    Ut(fidx) = -c*(U(fidx) -U(bidx))/dx;
    Ut(1) = Ut(end);
    
    YRK = [Ut; U];
    
    k1 = propGVD(t,YRK,options,A);
%     k1(1,:) = k1(N,:); k1(N+1,:) = k1(2*N,:);
    k2 = propGVD(t+dt/2,YRK+dt/2*k1,options,A);
%     k2(1,:) = k2(N,:); k2(N+1,:) = k2(2*N,:);
    k3 = propGVD(t+dt/2,YRK+dt/2*k2,options,A);
%     k3(1,:) = k3(N,:); k3(N+1,:) = k3(2*N,:);
    k4 = propGVD(t+dt,YRK+dt*k3,options,A);
%     k4(1,:) = k4(N,:); k4(N+1,:) = k1(2*N,:);
    
    YRK = YRK+dt/6*(k1+2*k2+2*k3+k4);
    Ut = YRK(1:N) ; U = YRK((N+1):end);
    
    solver.set_latest_solution(U);
    
%     U = solver.make_step(dummy,dummy,dummy,dt/2);
%     U = solver.set_bdry( U(end),'no');
%     
    if(mod(ctr,100) == 0)
        plot(x,abs(U),x,abs(Uorig));
        title([ 't = ' num2str(t) ] );
        getframe;
    end
    
    ctr =ctr+1;
    t=t+dt;
end





clear;clc; close all;

c = 1; n = 1; 
a = 0; b = 5;
beta = 1e5; N = 1000; 
E_0 = 1; initial_conditions = zeros(N,1);
x = linspace(a,b,N)';
dx = x(2)-x(1);
dt = 1*dx/(abs(c));

tau =0.01; %sec

aE = @(z,t) E_0.*sech((t-(z-(b+a)/2)/c)/tau);
aE_t =  @(z,t) -E_0/tau.*sech((t-(z-(b+a)/2)/c)/tau).*tanh((t-(z-(b+a)/2)/c)/tau);
% 
sig2 = ((b-a)/100)^2;
aE = @(z,t) E_0.*exp(-(t-z+(b+a)/2).^2/2/sig2);
aE_t =  @(z,t)-E_0*(t-z+(b+a)/2).*exp(-(t-z+(b+a)/2).^2/2/sig2) ;


U = aE(x,0); V = aE_t(x,0);


ctr = 1;
t = dt;

tend = 100*(2*(b-a)/abs(c));
U_t = 0 ;

fidx = 2:N; 
bidx = 1:N-1;
Unew = U; Vnew = V; 
while t<tend
    
    %backward differentiation
    Ux = (U(fidx)-U(bidx))/dx; 
    Vx = (V(fidx)-V(bidx))/dx;
    
    
    Vnew(fidx) = V(fidx) + dt*2*1i/beta*( n/c*V(fidx)+Ux) + dt^2/2*(-4*(n/beta/c)^2*V(fidx) - 4*n/beta^2/c*Ux +2*1i/beta*Vx);
    Unew(fidx) = U(fidx) + (dt^2/2)*2*1i/beta*(n/c*V(fidx)+Ux) + dt*V(fidx);
%     Vnew(1) =Vnew(end);
%     Unew(1) = Unew(end); 
    
    U = Unew; V = Vnew; 
    
    if(mod(ctr,1) == 0)
       
        plot(x,abs(U));
        title([ 't = ' num2str(t) ] );
        getframe;
        
    end
    
    
    
    ctr =ctr+1;
    t=t+dt;
end





clear;clc; close all;

c = 2 ; a = 0; b = 5;
beta = 0.01;

D = -1i*beta*c^3/4; 

N = 1000; dx = (b-a)/(N);
x = linspace(a-dx,b+dx,N)';

dx = x(2)-x(1); 
E_0 = 1; tau =0.01; %sec

aE = @(z,t) E_0.*sech((t-(z-(b+a)/2)/c)/tau);
aE_t =  @(z,t) -E_0/tau.*sech((t-(z-(b+a)/2)/c)/tau).*tanh((t-(z-(b+a)/2)/c)/tau);


Ep = aE(x,0);
Em = 0*Ep;


ctr = 1;

t = 0;

tend = 100*(2*(b-a)/abs(c));
options = '';
Et0 = 0 ; 
Etend = 0;

EPnew = Ep;
Energyt = 0;

sdt = 0.9*dx/abs(c);
while t<tend

   EPnew(2:end) = Ep(2:end)-c*(dt/dx)*(Ep(2:end)-Ep(1:end-1));
   EPnew(1) = EPnew(end); %periodic boundaries! 

%    plot(x,abs(EPnew));
%    getframe;

   Ep = EPnew; 
%    sum(abs(Ep).^2)*dx
   Energyt(ctr) = trapz(x,abs(Ep).^2);
   
   t= t+dt;
   ctr = ctr+1;
end





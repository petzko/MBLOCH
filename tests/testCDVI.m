clear;clc; close all;

c = +1; a = -2; b = 2; TR = (b-a)/2/c;
beta =1; N =2000; 
% x = (a+b)/2-(b-a)/2*cos(pi*(0:N-1)/(N-1))';
x = linspace(a,b,N)';
dx = x(2)-x(1); dt = dx/abs(c);
E_0 = 1;



tau =0.1; %sec
k = 1;
z0 = (b+a)/2
aE = @(z,t) E_0.*exp(-((t-(z-z0)/c)/tau).^2).*exp(1i*x*k);
aE_t =  @(z,t) -1/5.*2*(t-(z-z0)/c)/tau.* E_0.*exp(-((t-(z-z0)/c)/tau).^2).*exp(1i*x*k);;


U = aE(x,0); 
V = aE_t(x,0);


U(1) = 0;  U(end) = 0 ; V(1) = 0  ; V(end) = 0; 


ctr = 1;
t = 2*dt;
tend = 20*TR;

solver = RNFDSolver(N,dx,c/abs(c),abs(c), U);
dummy = 0*U;

while t<tend

    
    if(mod(ctr,100)==0)
        plot(x,abs(U));
        title(['t @ ' num2str(t)]);
        getframe;
    end
   
%     U = solver.make_step(dummy,dummy,dummy,dt/2);
%     U = solver.set_bdry(U(end),'no');
%     
    U = U + (dt/2+1i*dt^2/4/beta/c)*V;
    V = V.*exp(2*1i/beta/c*dt/2);
    solver.set_latest_solution(U);
    
    U = solver.make_step(dummy,dummy,dummy,dt);
    U = solver.set_bdry(U(end),'no');
    
    U = U + (dt/2+1i*dt^2/4/beta/c)*V;
    V = V.*exp(2*1i/beta/c*dt/2);
    
    solver.set_latest_solution(U);    
    ctr =ctr+1;
    t=t+dt;
end





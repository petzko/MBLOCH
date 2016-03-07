clear;close all;
% reflection coefficient

R = 0.75;  

nr_rt = 4;
a = -2;
b = 2; 
L = b-a; %% length of the domain
c = 0.6;  %% propagation velocity
T_R = 2*L/c ;  %% roundtrip time

tEnd = nr_rt*T_R;

N = 128;  %% nr of grid points 

t = 0;
x = linspace(a,b,N); 
dx = x(2) - x(1);
dt =0.01*dx/c; 
x = (a+b)/2 - (b-a)/2*cos(pi*(0:N-1)/(N-1));


Au = 3*exp(-(x-1).^2/(2*0.01))/sqrt(2*pi*0.01); 
Av = 0;
k = 1;
F_rk = (Au.*exp(1i*k*x))';
F_lax = F_rk;
F_spec = F_rk; 

B_rk = (Av.*exp(1i*k*x))';
B_lax = B_rk;
B_spec = B_rk; 


Fspec_solver = SpectralMSolver(N,1,0,L,3,c,[],F_spec);
Flax_solver =  LaxSolver(N,dx,1,c, F_lax);

Bspec_solver = SpectralMSolver(N,-1,0,L,3,c,[],B_spec);
Blax_solver =  LaxSolver(N,dx,-1,c, B_lax);

idx = round(N/2);
ampl(1) = F_rk(idx).*conj(F_rk(idx)) + B_rk(idx).*conj(B_rk(idx)); 
ctr =  2; 
skip = 10;
iter_ctr = 1;
dummy = zeros(N,1);

%chebishev differentiation matrix! 
Dmtx = 2/L*cheb(N-1);

bp = zeros(N,1); bp(1) = 1; 
bn = zeros(N,1); bn(end) = 1; 

rhsp = @(f) -c*Dmtx*f;
rhsn = @(f) +c*Dmtx*f;

while(t < tEnd)
 
        
        % first kind of bdry conditions
        
          k1 = rhsp(F_rk);
%           k1(1) = R*k1(end); 
          
          k2 = rhsp(F_rk + dt/2*k1);
%           k2(1) = R*k2(end); 
 
          
          k3 = rhsp(F_rk + dt/2*k2);
%           k3(1) = R*k3(end);

          k4 = rhsp(F_rk + dt*k3);
%           k4(1) = R*k4(end);

          F_rk = F_rk + (dt/6)*(k1+2*k2+2*k3+k4);
          F_rk(1) = R*F_rk(end);
          
          F_lax = Flax_solver.make_step(dummy,dummy,dummy,dt);
%           B_lax = Blax_solver.make_step(dummy,dummy,dummy,dt);
          
          F_lax = Flax_solver.set_bdry(R*F_lax(end),'no');
%           B_lax = Blax_solver.set_bdry('no',F_lax(end));

          F_spec = Fspec_solver.make_step(dummy,dt);
%           B_spec = Bspec_solver.make_step(dummy,dt);
          
          F_spec = Fspec_solver.set_bdry(R*F_spec(end),'no');
%           F_spec = Fspec_solver.set_bdry(R*B_spec(1),'no');
%           B_spec = Bspec_solver.set_bdry('no',R*F_spec(end));
          
          E_rk = F_rk.*conj(F_rk) + B_rk.*conj(B_rk);
          E_lax = F_lax.*conj(F_lax) + B_lax.*conj(B_lax);
          E_spec = F_spec.*conj(F_spec) + B_spec.*conj(B_spec);
          
      
          
        if(mod(iter_ctr,skip) == 0 )
            ampl(ctr) = E_rk(idx);
            ctr = ctr+ 1;
        end
      if(mod(iter_ctr,100) == 0 )
        plot(x,E_rk,x,E_lax,x,E_spec);
        axis([a b 0 200])
        legend('rk', 'lax', 'spec')
        getframe;
      end
        iter_ctr = iter_ctr + 1;
        
    t = t + dt;
end



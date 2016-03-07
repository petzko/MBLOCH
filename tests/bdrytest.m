R = 1;
nr_rt = 4;
a = -1;
b = 4;

c = 1;
T_R = 2*(b-a)/c ;

tEnd = nr_rt*T_R;

N = 200;

t = 0;
x = linspace(a,b,N);
x_spec = (a+b)/2-(b-a)/2*cos(pi*(0:(N-1))/(N-1));
dx = x(2) - x(1);
dt = 0.01*dx/c;

dx_rnfd = c*dt;
N_rnfd = round((b-a)/dx_rnfd);  
x_rnfd =  linspace(a,b,N_rnfd); % -cos(pi*(0:(N_rnfd-1))/(N_rnfd-1))'; 
dx_rnfd = x_rnfd(2) - x_rnfd(1);
dt_rnfd = dx_rnfd/c;


A_0 = 3*exp(-(x_spec-0.75*(b-a)).^2/(2*0.01))/sqrt(2*pi*0.01);
A_0_rnfd = 3*exp(-(x_rnfd-0.75*(b-a)).^2/(2*0.01))/sqrt(2*pi*0.01);
A_1 = 0;
k = 1;
u_lax = (A_0.*exp(1i*k*x_spec))';
v_lax = (A_1.*exp(1i*k*x_spec))';
u_rnfd = (A_0_rnfd.*exp(1i*k*x_rnfd))';
v_rnfd = (A_1.*exp(1i*k*x_rnfd))';
u_horder = u_lax;
v_horder = v_lax;

a1 = c*dt/(2*dx);
b1 = 2*a1*a1;


idx = round(N/2);
ampl(1) = u_lax(idx).*conj(u_lax(idx)) + v_lax(idx).*conj(v_lax(idx));
ctr =  2;
skip = 10;
iter_ctr = 1;
dummy = zeros(N,1);

F_lax = LaxSolver(N ,dx,1,c, u_lax);
B_lax = LaxSolver(N ,dx,-1,c, v_lax);

F_rnfd = RNFDSolver(N_rnfd ,dx_rnfd,1,c, u_rnfd);
B_rnfd = RNFDSolver(N_rnfd ,dx_rnfd,-1,c, v_rnfd);



F_HO = SpectralMSolver(N,1,a,b,3,c,[],u_horder);
B_HO = SpectralMSolver(N,-1,a,b,3,c,[],v_horder);



dummy = zeros(N,1); dummy_rnfd = zeros(N_rnfd,1);

while(t < tEnd)
    
    
    
    %%% lax
    u_lax = F_lax.make_step(dummy,dummy,dummy,dt);
    v_lax = B_lax.make_step(dummy,dummy,dummy,dt);
    
    u_left = R*v_lax(1);
    v_right= R*u_lax(end);
    
    u_lax = F_lax.set_bdry( u_left ,'no' );
    v_lax = B_lax.set_bdry('no',v_right);
    
    %%% rnfd
    u_rnfd = F_rnfd.make_step(dummy_rnfd,dummy_rnfd,dummy_rnfd,dt_rnfd);
    v_rnfd = B_rnfd.make_step(dummy_rnfd,dummy_rnfd,dummy_rnfd,dt_rnfd);
    
    u_left = R*v_rnfd(1); v_right= R*u_rnfd(end);
    
    u_rnfd = F_rnfd.set_bdry( u_left ,'no' );
    v_rnfd = B_rnfd.set_bdry('no',v_right);
    
    %%% spectral
    u_horder = F_HO.make_step(dummy,dt);
    v_horder = B_HO.make_step(dummy,dt);
    
    u_left = R*v_horder(1); 
    v_right = R*u_horder(end);
    
    u_horder = F_HO.set_bdry( u_left ,'no' );
    v_horder = B_HO.set_bdry('no',v_right);
    
    E_ho = u_horder.*conj(u_horder) + v_horder.*conj(v_horder);
    E_lax = u_lax.*conj(u_lax) + v_lax.*conj(v_lax);
    E_rnfd = u_rnfd.*conj(u_rnfd) + v_rnfd.*conj(v_rnfd);

    if(mod(iter_ctr,skip) == 0 )
        ampl(ctr) = E_lax(idx);
        ctr = ctr+ 1;
    end
    
    if(mod(iter_ctr,10) == 0)
        plot(x,E_lax,x_rnfd,E_rnfd,x,E_ho);
        legend('lax','rnfd','ho');
        axis([a b 0 200])
        getframe;
    end
    iter_ctr = iter_ctr + 1;
    
    t = t + dt;
end



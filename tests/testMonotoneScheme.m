L = 1;  N = 256;
x = linspace(0,L,N)';
% x = L/2-L/2*cos(pi*(0:N-1)/(N-1))';
dx = x(2) - x(1); 
c = 0.25; tend = 1.5; 
dt = 0.7*dx/c;

E_f = zeros(N,1); E_b = zeros(N,1); 
init = round(0.25/dx);
fin = round(0.5/dx); 
xsmall = linspace(0,(fin-init)*dx,fin-init+1)';
k = 2*2*pi/xsmall(end);
E_f(init:fin) = 0.5*sin((xsmall-init*dx)*k);
E_finit = E_f;

forward_wave = LaxSolver(N,dx,1,c, E_f);
backward_wave = LaxSolver(N,dx,-1,c, zeros(N,1));

% forward_wave = SpectralMSolver(N,1,0,L,3,c,[],E_f);
% backward_wave = SpectralMSolver(N,-1,0,L,3,c,[],E_b);

t = dt;
F = zeros(N,1); F_t = F; k = F;
iter = 0;

while t < tend
   
 
    E_f = forward_wave.make_step(F,F_t,k,dt);
    E_b = backward_wave.make_step(F,F_t,k,dt);
%     E_f = forward_wave.make_step(F,dt);
%     E_b = backward_wave.make_step(F,dt);

    E_f = forward_wave.set_bdry(E_b(1),'no');
    E_b = backward_wave.set_bdry('no',E_f(end));     
    
    t = t+dt; 
    iter = iter + 1; 
end
        figure;
       ax =  plot(x,E_finit,'-b',x,E_f,'--r');
        title(['t = ' num2str(t)])
        axis([0 L -0.7 0.7]);
    %     set(ax(1),'Ylim',[-0.1 1.5],'Line','--r'); 
    %     set(ax(2),'Ylim',[-0.1 1.5]); 
        getframe;
    
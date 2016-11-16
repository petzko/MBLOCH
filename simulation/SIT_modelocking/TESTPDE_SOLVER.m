% length of the cavity  
L = 3;
% propagation speed 
c = .1; 
% number of points for spatial discretization 
N = 3000;
x = linspace(0,L,N)';
%grid spacing
dx = x(2)-x(1);
% time step -> Always choose dt = dx/c;  
dt = dx/c;
% end of simulation time 
tEnd = 100; 

x0 = L/2
sigma = 0.1;
U = exp(-(x-x0).^2/2/sigma^2); 

% initialize PDE solver (RNFDSolver) with IC = U, velocity = c and 
% propagating in the positive x direction +1 
solver = RNFDSolver(N,dx,+1,c,U)

t = dt; 
dummy = zeros(N,1);
loss = -1e-2*ones(N,1); 

iter_ctr = 1;
while t<tEnd
    U = solver.make_step(dummy,dummy,loss,dt);
%     U = solver.set_bdry(.3,'no'); % stupid boundary conditions
    U = solver.set_bdry(U(end),'no'); % periodic boundary conditions
    
    if mod(iter_ctr,50) == 0
        plot(x,abs(U))
        ylim([0,1])
        getframe;
    end

    t = t+dt;
    iter_ctr = iter_ctr+1;
        
end





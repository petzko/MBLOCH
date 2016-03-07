%%%% examine the convergence properties of all currently implemented
%%%% hyperbolic equation solvers! 

c = 1; %forward wave...
U = @(x,t) sin(2*pi*(x/c-t)) + 5;
tend =10; 
L =3; 

N = 20;
max_attempts = 2;
max_error = zeros(max_attempts,1);

for k = 1:max_attempts
    N = N*2;
    x = linspace(0,L,N);
    dx = x(2) - x(1);
    U_num = U(x,0)';
    
%     dt = 0.5*dx/c;
%     solver = LaxSolver(N,dx,1,c,U_num);
    dt = dx/c;
    solver = RNFDSolver(N,dx,1,c,U_num);
%     dt = 0.5*dx/c;
%     solver = HighOrderSolver(N,dx,dt,1,c,U_num);

    errors = 0; 
    ctr = 0;
    t = 0;
    display(['attempt ' num2str(k) '; N = ' num2str(N)]);
    while(t < tend) 
        ctr = ctr + 1;
        error_x= abs(U(x,t)' - U_num);
        errors(ctr) = max(error_x);
        plot(x,U(x,t),x,U_num);
        title(['time = ' num2str(t)])
        axis([0 L 4 6]);
        getframe;
%         U_num = solver.make_step(0,zeros(N,1),zeros(N-2,1),zeros(N-2,1),dt);
        U_num = solver.make_step(zeros(N,1),zeros(N,1),0,dt);
%         U_num = solver.make_step(zeros(N,1),dt);
        U_num = solver.set_bdry(U(x(1),t+dt),'no');
        t =t+dt;
    end
    display(max(errors));
    max_error(k,1) = max(errors);
end
% 
plot(1:max_attempts,max_error)
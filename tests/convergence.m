clear;clc;
%analytic solution to our problem...
%%%%%%%%%% Spatial error order measurement


w = 4;
v = -1;
k = w/v;
a = 0;
b = 3;
%analytic solution to our problem...
constant = 1.5;
u_analitic = @(x,t) (x*t.*exp(1i*(k*x-w*t)) + constant);
f = @(x,t) (x + v*t).*exp(1i*(k*x-w*t));
f_t = @(x,t) v*exp(1i*(k*x-w*t)) - 1i*w*(x+v*t).*exp(1i*(k*x-w*t));
f_x = @(x,t) exp(1i*(k*x-w*t)) + 1i*k*(x+v*t).*exp(1i*(k*x-w*t));
% u_analitic = @(x,t) (sin(w*(x-t/v)));
% f = @(x,t) 0*x;
% f_t = @(x,t) 0*x;
% f_x = @(x,t) 0*x;
plot_ctr = 2500;
tend = 10;

%% dx- convergence general setup


%%%%%%%%%%%%%%%
orderend = 3;
N_init = 64;
sizes(1) = N_init;
for i =2:orderend
    sizes(i) = 2*sizes(i-1);
end
sizes = sizes';
% simulation

%arrays to store the global error or euler and lax method for each dx configuration
dx_error_samples = zeros(orderend,1);
for order = 1:orderend
    %% for loop dx - convergence
    nr_steps = 3;
    display(['space error measurement : order = ' num2str(order)]);
    N = sizes(order); % nr of grid points
    display(['nr of grid points ' num2str(N)]);
    x = (a+b)/2-(b-a)/2*cos(pi*(0:N-1)/(N-1))';
%     x = linspace(a,b,N);
    dx = x(2)-x(1);
    dt = 1.2*dx/abs(v);
    %     DX = x(2:end-1) - x(1:end-2);
    
    [ M, N ] = size(x);
    N = max(M,N);
    initial_conditions = zeros(N,nr_steps-1);
    
    %         %prepare the initial conditions...
    for k = 1:nr_steps-1
        initial_conditions(:,k) = f(x,(k-1)*dt);
    end
    %initialize the analytic solution
    y_analytic = u_analitic(x,(nr_steps-1)*dt);
    u_num = u_analitic(x,(nr_steps-1)*dt);
%     solver = RNFDSolver(N,dx,v/abs(v),abs(v), u_num);
%     solver = LaxSolver(N,dx,v/abs(v),abs(v), u_num);
    %       solver = HighOrderSolverImplicit(N,dx,dt,direction,velocity, U_0)
    %       solver = HighOrderSolver(N,dx,dt,direction,velocity, U_0)
    solver = SpectralMSolver(N,v/abs(v),a,b,nr_steps,abs(v),initial_conditions,u_num);
    
    t = (nr_steps-1)*dt;
    ctr = 1;
    num_err = 0;
    
    while(t <= tend)
        
        % main loop dx convergence...
                if(mod(ctr,10) == 0)
                    plot(x,real(u_num),'g',x,real(y_analytic),'--b');
                %    axis([a b -1 1])
                %   legend('lax-wendroff ','euler','analytic sol');
                    title(['accuracy comparison: @ t = ' num2str(t) '; error = ' num2str(max(num_err))])
                    getframe;
                end
        %
        num_err(ctr) = 0;
        %analytic solution
        y_analytic = u_analitic(x,t+dt);
% %         u_num = solver.make_step(f(x,t),f_t(x,t),0,dt);
        u_num = solver.make_step(f(x,t), dt);
%                 u_num = solver.set_bdry(y_analytic(1),y_analytic(end));
        if(v < 0)
            u_num = solver.set_bdry('no',y_analytic(end));
        else
            u_num = solver.set_bdry(y_analytic(1),'no');
        end
        %error estimate
        num_err(ctr) =(max(u_num - y_analytic));
        t = t+dt;
        ctr = ctr+1;
        
    end
% % %     store the maximum errors...
        dx_error_samples(order,:) = [max(num_err)];
    
end
%     figure;
%     plot(sizes,dx_error_samples(:,1),'r');
%     xlabel('N (nr. grid points)' )
%     ylabel('|\epsilon| - global error');
%     legend(['\epsilon_{num}'] );
%     title([' RNFD Convergence with \omega = ' num2str(w)]);
display(abs(dx_error_samples))


clear;clc;
%analytic solution to our problem...
%%%%%%%%%% Spatial error order measurement

w = 4;
v = +1;

k = w/v;
a = 0;
b = 3;
%analytic solution to our problem...
constant = 1.5;
u_analitic = @(x,t) (x*t.*exp(1i*(k*x-w*t)) + constant);
f = @(x,t) (x + v*t).*exp(1i*(k*x-w*t));
f_t = @(x,t) v*exp(1i*(k*x-w*t)) - 1i*w*(x+v*t).*exp(1i*(k*x-w*t));
f_x = @(x,t) exp(1i*(k*x-w*t)) + 1i*k*(x+v*t).*exp(1i*(k*x-w*t));
plot_ctr = 1;
tend = 5;
x_measure = 1 ; t_measure = -1;

%% dx- convergence general setup
if(x_measure > 0)
    
    %%%%%%%%%%%%%%%
    orderend = 8;
    N_init = 64;
    sizes(1) = N_init;
    for i =2:orderend
        sizes(i) = 2*sizes(i-1);
    end
    N_end = sizes(end);
    sizes = sizes';
    % simulation
    
    %arrays to store the global error for euler and lax method for each dx configuration
    dx_error_samples = zeros(orderend,1);
    
    for order = 1:orderend
        
        display(['global error measurement : order = ' num2str(order)]);
        N = sizes(order); % nr of grid points
        x = linspace(a,b,N)';
        dx = x(2)-x(1);
        dt = dx/abs(v);
        
        y_analytic = u_analitic(x,0);
        u_euler = u_analitic(x,0);
        %initialize lax and some parameters... x
        u_num = u_analitic(x,0);
        solver = RNFDSolver(N,dx,v/abs(v),abs(v), u_num);
        t = dt;
        ctr = 1;
        num_err = 0;
        
        while(t <= tend)
            if(mod(ctr,1000) == 0)
                plot(x,real(u_num),'g',x,real(y_analytic),'--b');
                title(['Numerical and analytic solution @ t = ' num2str(t)])
                getframe;
            end
            
            num_err(ctr) = 0;
            %analytic solution
            
            y_analytic = u_analitic(x,t);
            u_num = solver.make_step(f(x,t-dt),f_t(x,t-dt),zeros(N,1),dt);
            
            if(v < 0)
                u_num = solver.set_bdry('no',y_analytic(end));
            else
                u_num = solver.set_bdry(y_analytic(1),'no');
            end
            
            %error estimate
            num_err(ctr) =norm(u_num - y_analytic);
            t = t+dt;
            ctr = ctr+1;
        end
        %store the maximum errors...
        dx_error_samples(order,:) = [max(num_err)];
        
    end
end

plot(2:orderend,dx_error_samples(1:end-1)./dx_error_samples(2:end))


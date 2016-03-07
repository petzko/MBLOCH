%analytic solution to our problem...
x_measure = +1 ; t_measure = -1
%%%%%%%%%% Spatial error order measurement
if(x_measure > 0)
    %%%%%%%%%%%^
    w = 5;
    v = 1;
    k = w/v;
    %analytic solution to our problem...
    u_analitic = @(x,t) (x*t.*exp(1i*(k*x-w*t)));
    f = @(x,t) (x + v*t).*exp(1i*(k*x-w*t));
    f_t = @(x,t) v*exp(1i*(k*x-w*t)) - 1i*w*(x+v*t).*exp(1i*(k*x-w*t))
    f_x = @(x,t) exp(1i*(k*x-w*t)) + 1i*k*(x+v*t).*exp(1i*(k*x-w*t));
    xend = 10;
    tend = 0.5;
    %%%%%%%%%%%%%%%
    orderend = 50;
    increment = 20;
    dx_samples = [increment:increment:orderend*increment];
    dt = 0.1*xend/(v*orderend*increment-1);
    dt_error_samples = zeros(orderend,2); %% arrays to store the global error for euler and lax method for each dt configuration
    %measure the error when we gradually decrease the grid size...
    for order = 1:orderend;
        clc;
        display(['space error measurement : order = ' num2str(order)]);
        N = order*increment; % nr of grid points will be 50,100,150.. 1000
        x= linspace(0,xend,N)';
        dx = x(2)-x(1);
        %boundary condition at left boundary
        g_l = 0;
        %boundary condition at right boundary
        g_r = @(t) u_analitic(xend,t);
        %initial conditions...
        u_numeric = zeros(size(x));
        [ M, N ] = size(x);
        N = max(M,N);
        y_analytic = zeros(N,1);
        u_euler = zeros(N,1);
        %initialize lax and some parameters... x
        u_lax = zeros(N,1);
        y_x = zeros(N-2,1);
        y_xx = zeros(N-2,1);
        a1 = v*dt/(2*dx);
        b1 = 2*a1*a1;
        c1 = dt;
        d1 = (dt*dt)/2;
        x_short = x(2:N-1);
        t = dt;
        ctr = 1;
        while(t <= tend)
            err_euler(ctr) = 0;
            err_lax(ctr) = 0;
            %analytic solution
            y_analytic = u_analitic(x,t);
            %numerical solution using euler
            tmp = f(x,t);
            u_euler(2:N) = u_euler(2:N) - v*dt/dx*(u_euler(2:N) - u_euler(1:N-1)) + dt*tmp(2:N);
            u_euler(1) = y_analytic(1);
            %numerical solution using lax method...
            y_x = u_lax(3:N) - u_lax(1:N-2);
            y_xx = u_lax(3:N) - 2*u_lax(2:N-1) + u_lax(1:N-2);
            u_lax(2:N-1,1) = u_lax(2:N-1,1) - a1*y_x + b1*y_xx +c1*f(x_short,t) + d1*(f_t(x_short,t) - v*f_x(x_short,t));
            u_lax(1) = y_analytic(1);
            u_lax(end) = y_analytic(end);
            plot(x,real(u_lax),'--r',x,real(y_analytic),'g',x,real(u_euler),'b');
            getframe
            %error estimate
            err_euler(ctr) = max(abs(u_euler-y_analytic));
            err_lax(ctr) = max(abs(u_lax - y_analytic));
            t = t+dt;
            ctr = ctr+1;
        end
        %store the maximum errors...
        dx_error_samples(order,:) = [max(err_euler) max(err_lax)];
    end
    C = 10*20000;
    f = C./dx_samples.^2 ;
    plot(dx_samples,dx_error_samples(:,1),'r',dx_samples,dx_error_samples(:,2),'g',dx_samples,f);
    axis([0 dx_samples(end) 0 50])
    xlabel('N (nr. grid points)' )
    ylabel('|\epsilon| - global error');
    legend(['\epsilon_{euler}'],['\epsilon_{lax}'], ['C/N^2 , with C = ' num2str(C) ] );
    title(['Lax-Wendroff scheme convergence dt fixed @ ' num2str(dt) ' ; dx varying']);
end
if(t_measure > 0)
    w = 20;
    v = 1;
    k = w/v;
    %analytic solution to our problem...
    u_analitic = @(x,t) (x*t.*exp(1i*(k*x-w*t)));
    f = @(x,t) (x + v*t).*exp(1i*(k*x-w*t));
    f_t = @(x,t) v*exp(1i*(k*x-w*t)) - 1i*w*(x+v*t).*exp(1i*(k*x-w*t))
    f_x = @(x,t) exp(1i*(k*x-w*t)) + 1i*k*(x+v*t).*exp(1i*(k*x-w*t));
    xend = 10;
    tend = 0.5;
    %%%% MEASURE THE ERROR IN THE TIME COORDINATE
    orderend = 50;
    dt_init = 10^-3; % initial timestep is 0.009...
    dt_samples = linspace(dt_init,dt_init/orderend,orderend);
    dx = 1.1*dt_init; %%initial dx is 0.0099
    x = 0:dx:xend; x= x'; %transpose x
    dt_error_samples = zeros(orderend,2); %% arrays to store the global error for euler and lax method for each dt configuration
    %measure the error when we gradually decrease timestep..
    for order = 1:orderend;
        clc;
        display(['time error measurement : order = ' num2str(order)]);
        dt = dt_samples(order);
        %boundary condition at left boundary
        g_l = 0;
        g_r = @(t) u_analitic(xend,t);
        %boundary condition at right boundary
        %initial conditions...
        u_numeric = zeros(size(x));
        [ M, N ] = size(x);
        N = max(M,N);
        y_analytic = zeros(N,1);
        u_euler = zeros(N,1);
        %initialize lax and some parameters... x
        u_lax = zeros(N,1);
        y_x = zeros(N-2,1);
        y_xx = zeros(N-2,1);
        a1 = v*dt/(2*dx);
        b1 = 2*a1*a1;
        c1 = dt;
        d1 = (dt*dt)/2;
        x_short = x(2:N-1);
        t = 0;
        ctr = 1;
        while(t <= tend)
            err_euler(ctr) = 0;
            err_lax(ctr) = 0;
            %analytic solution
            y_analytic = u_analitic(x,t);
            %numerical solution using euler
            tmp = f(x,t);
            u_euler(2:N) = u_euler(2:N) - v*dt/dx*(u_euler(2:N) - u_euler(1:N-1)) + dt*tmp(2:N);
            u_euler(1) = y_analytic(1);
            %numerical solution using lax method...
            y_x = u_lax(3:N) - u_lax(1:N-2);
            y_xx = u_lax(3:N) - 2*u_lax(2:N-1) + u_lax(1:N-2);
            u_lax(2:N-1,1) = u_lax(2:N-1,1) - a1*y_x + b1*y_xx +c1*f(x_short,t) + d1*( f_t(x_short,t) - v*f_x(x_short,t) );
            u_lax(1) = y_analytic(1);
            u_lax(end) = y_analytic(end);
            %error estimate
            err_euler(ctr) = max(abs(u_euler-y_analytic));
            err_lax(ctr) = max(abs(u_lax - y_analytic));
            t = t+dt;
            ctr = ctr+1;
        end
        %store the maximum errors...
        dt_error_samples(order,:) = [max(err_euler) max(err_lax)];
    end
    C = -10000;
    B = 150;
    f = C*dt_samples.^2 + B*dt_samples ;
    figure; % plot in a new figure
    plot(dt_samples,dt_error_samples(:,1),'r',dt_samples,dt_error_samples(:,2),'g',dt_samples,f);
    % axis([0 dt_samples(end) 0 50]);
    xlabel('dt (timestep)');
    ylabel('|\epsilon| - global error');
    legend(['\epsilon_{euler}'],['\epsilon_{lax}'], ['c*dt^2 + b*dt ,(c = ' num2str(C) ' ; b = ' num2str(B) ')']);
    title(['Lax-Wendroff scheme convergence dx fixed @ ' num2str(dx) ' ; dt varying']);
end
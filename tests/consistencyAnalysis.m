clear;clc;
%analytic solution to our problem...
%%%%%%%%%% Spatial error order measurement

w = 10;
c = -1;
k = w/c;
a = -1;
b = 1;
%analytic solution to our problem...
u_analitic = @(x,t) (x*t.*exp(1i*(k*x-w*t)));
f = @(x,t) (x + c*t).*exp(1i*(k*x-w*t));
f_t = @(x,t) c*exp(1i*(k*x-w*t)) - 1i*w*(x+c*t).*exp(1i*(k*x-w*t));
f_x = @(x,t) exp(1i*(k*x-w*t)) + 1i*k*(x+c*t).*exp(1i*(k*x-w*t));
plot_ctr = 1000;
tend = 0.5;
% which convergence criteria to measure
x_measure = 1 ; t_measure = -1;
    


%% dx- convergence general setup
if(x_measure > 0)
    
    %%%%%%%%%%%%%%%
    orderend = 5;
    N_init = 64;
    sizes(1) = N_init; 
    for i =2:orderend
        sizes(i) = 2*sizes(i-1);
    end
    N_end = sizes(end);
    sizes = sizes';
    % simulation 
    
    %arrays to store the global error for euler and lax method for each dx configuration
    dx_error_samples = zeros(orderend,3); 
    
%   dt = 0.5*((b-a)/sizes(end))/abs(v);
    for order = 1:orderend;
        %% for loop dx - convergence
        clc;
        display(['space error measurement : order = ' num2str(order)]);        
        N = sizes(order); % nr of grid points
        display(['nr of grid points ' num2str(N)]);
        x = [linspace(a,b,N); linspace(a,b,N); (a+b)/2-((b-a)/2)*cos(pi*(0:N-1)/(N-1))];
        dx = [x(1,2)-x(1,1);x(2,2) - x(2,1); x(3,2)-x(3,1)];
        dt = [dx(1);0.5*dx(2);1.1*dx(3)]/abs(c);
    
        %initialize lax and some parameters... x
        u_num = [u_analitic(x(1,:),0);u_analitic(x(2,:),0);u_analitic(x(3,:),0)];
        rnfd_solver = RNFDSolver(N,dx(1),c/abs(c),abs(c), u_num(1,:)');
        lax_solver  = LaxSolver(N,dx(2),c/abs(c),abs(c), u_num(2,:)');
        spectral_solver = SpectralMSolver(N,c/abs(c),a,b,3,abs(c),u_num(3,:)');
        t = [dt(1);dt(2);dt(3)];
        ctr = [1;1;1];
        errors_rnfd = 0;
        errors_lax = 0;
        errors_spectral =0;
        %% main loop dx convergence RNFD... 
        while(t(1) <= tend)
%             if(mod(ctr,plot_ctr) == 0)
%                 plot(x(1,:),real(u_num(1,:)),'g',x(1,:),real(u_analitic(x(1,:),t(1,:)-dt(1))),'--b');
%                 axis([a b -1 1])
%                 title(['analytic vs numerical solution @ t = ' num2str(t(1))])
%                 getframe;
%             end
            
            % new analytic solution
            y_analytic = u_analitic(x(1,:),t(1));
            u_num(1,:) = rnfd_solver.make_step(f(x(1,:),t(1)-dt(1))',f_t(x(1,:),t(1)-dt(1))',0,dt(1));
           
            if c > 0
                u_num(1,:) = rnfd_solver.set_bdry(y_analytic(1),'no');
            else
                u_num(1,:) = rnfd_solver.set_bdry('no',y_analytic(end));
            end
                %error estimate
            errors_rnfd(ctr(1)) = max(abs((u_num(1,:) - y_analytic)./abs(y_analytic)));
            t(1) = t(1)+dt(1);
            ctr(1) = ctr(1)+1;
        end
        

        %% main loop dx convergence LAX... 
        while(t(2) <= tend)
%             if(mod(ctr,plot_ctr) == 0)
%                 plot(x(2,:),real(u_num(2,:)),'g',x(2,:),real(u_analitic(x(2,:),t(2,:)-dt(2))),'--b');
%                 axis([a b -1 1])
%                 title(['analytic vs numerical solution @ t = ' num2str(t(2))])
%                 getframe;
%             end
            
            % new analytic solution
            y_analytic = u_analitic(x(2,:),t(2));
            u_num(2,:) = lax_solver.make_step(f(x(2,:),t(2)-dt(2))',f_t(x(2,:),t(2)-dt(2))',0,dt(2));
            u_num(2,:) = lax_solver.set_bdry(y_analytic(1),y_analytic(end));

            %error estimate
            errors_lax(ctr(2)) = max(abs((u_num(2,:) - y_analytic)./abs(y_analytic)));
            t(2) = t(2)+dt(2);
            ctr(2) = ctr(2)+1;
        end
        
        %% main loop dx convergence SPECTRAL
        while(t(3) <= tend)
%             if(mod(ctr,plot_ctr) == 0)
%                 plot(x(3,:),real(u_num(3,:)),'g',x(3,:),real(u_analitic(x(3,:),t(3,:)-dt(3))),'--b');
%                 axis([a b -1 1])
%                 title(['analytic vs numerical solution @ t = ' num2str(t(3))])
%                 getframe;
%             end
            
            % new analytic solution
            y_analytic = u_analitic(x(3,:),t(3));
            u_num(3,:) = spectral_solver.make_step(f(x(3,:),t(3)-dt(3))',dt(3));
            if c > 0
                u_num(3,:) = spectral_solver.set_bdry(y_analytic(1),'no');
            else
                u_num(3,:) = spectral_solver.set_bdry('no',y_analytic(end));
            end

            %error estimate
            errors_spectral(ctr(3)) = max(abs((u_num(3,:) - y_analytic)./abs(y_analytic)));
            t(3) = t(3)+dt(3);
            ctr(3) = ctr(3)+1;
        end
        
        
        %store the maximum errors...
        dx_error_samples(order,:) = [max(errors_rnfd) max(errors_lax) max(errors_spectral)];
        
    end
    
    plot(sizes,dx_error_samples);
    xlabel('N (nr. grid points)' )
     ylabel('|\epsilon| - global error');
    legend(['\epsilon_{num}'] );
    title(['Numerical scheme convergence analysis']);
end

%% dt-convergence measurement region
if(t_measure > 0)
    clc;
    %%%% MEASURE THE ERROR IN THE TIME COORDINATE
    orderend = 10;
    dt_init = 10^-1; % initial timestep is 0.009...
    dt_samples(1) = dt_init;
    for i = 2:orderend 
         dt_samples(i) = dt_samples(i-1)/2;
    end
    dt_samples = dt_samples';
    
    
    dx = 1.1*dt_init; %%initial dx is 0.0099
    x = 0:dx:xend; x= x'; %transpose x
    dt_error_samples = zeros(orderend,2); %% arrays to store the global error for euler and lax method for each dt configuration
    
    
    %measure the error when we gradually decrease timestep..
    for order = 1:orderend;
        %% dt - convergence 1 case 
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
        u_num = zeros(N,1);
        solver = LaxSolver(N,dx,1,0,c, u_euler);

        
        x_short = x(2:N-1);
        t = 0;
        ctr = 1;
        %% end of dt-convergence 1 case setup 
        while(t <= tend)
            %% dt convergence main loop
            err_euler(ctr) = 0;
            err_lax(ctr) = 0;
            %analytic solution
            y_analytic = u_analitic(x,t);
            %numerical solution using euler
            tmp = f(x,t);
            u_euler(2:N) = u_euler(2:N) - c*dt/dx*(u_euler(2:N) - u_euler(1:N-1)) + dt*tmp(2:N);
            u_euler(1) = y_analytic(1);
            %numerical solution using lax method...
            df_dt = (f(x_short,t) - f(x_short,t-dt))/dt;
            df_dx = (f(x_short,t) - f(x(1:N-2),t)/dx);
            solver.make_step(f(x,t),df_dt,df_dx,dt)
            u_num = solver.set_bdry(y_analytic(1));  
            if(mod(ctr,plot_ctr) == 0)
                plot(x,real(u_num),'--r',x,real(u_euler),'--g',x,real(y_analytic),'b');
                legend('lax-wendroff ','euler_sol','analytic sol');
                title(['accuracy comparison: @ t = ' num2str(t)])
                getframe;
            end
            %error estimate
            err_euler(ctr) = max(abs(u_euler-y_analytic));
            err_lax(ctr) = max(abs(u_num - y_analytic));
            t = t+dt;
            ctr = ctr+1;
        end
        %store the maximum errors...
        dt_error_samples(order,:) = [max(err_euler) max(err_lax)];
        
    end
  
    figure; % plot in a new figure
    plot(dt_samples,dt_error_samples(:,1),'r',dt_samples,dt_error_samples(:,2),'g');
    %     axis([0 dt_samples(end) 0 50]);
    xlabel('dt (timestep)');
    ylabel('|\epsilon| - global error');
    legend(['\epsilon_{euler}'],['\epsilon_{lax}']);
    title(['Lax-Wendroff scheme convergence dx fixed @ ' num2str(dx) ' ; dt varying']);
    
end


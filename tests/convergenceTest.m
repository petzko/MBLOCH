
clear;clc;
%analytic solution to our problem...
%%%%%%%%%% Spatial error order measurement

w = 4;
clear;clc;
%analytic solution to our problem...
%%%%%%%%%% Spatial error order measurement
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
% u_analitic = @(x,t) (sin(2*pi*(x-t)));
% f = @(x,t) 0*x;
% f_t = @(x,t) 0*x;
% f_x = @(x,t) 0*x;
plot_ctr = 1;
tend = 1;
% which convergence criteria to measure
x_measure = 1 ; t_measure = -1;

%% dx- convergence general setup
if(x_measure > 0)
    
    %%%%%%%%%%%%%%%
    orderend = 1;

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
    
    dt = 1.1*((b-a)/sizes(end))/abs(v);
    for order = 1:orderend
        %% for loop dx - convergence
        
        display(['space error measurement : order = ' num2str(order)]);
        N = sizes(order); % nr of grid points
        display(['nr of grid points ' num2str(N)]);
        x = (a+b)/2-(b-a)/2*cos(pi*(0:N-1)/(N-1))';
%         x = linspace(a,b,N);
        dx = x(2)-x(1);
        
        %         DX = x(2:end-1) - x(1:end-2);

    % simulation 
    
    %arrays to store the global error for euler and lax method for each dx configuration
    dx_error_samples = zeros(orderend,1); 
    compute_times = zeros(orderend,1);
    dt = 1.1*((b-a)/sizes(end))/abs(v);
    for order = 1:orderend
        %% for loop dx - convergence

        display(['space error measurement : order = ' num2str(order)]);        
        N = sizes(order); % nr of grid points
        display(['nr of grid points ' num2str(N)]);
        x = (a+b)/2-(b-a)/2*cos(pi*(0:N-1)/(N-1))'; 
%         x = linspace(a,b,N);
        dx = x(2)-x(1);
        DX = x(2:end-1) - x(1:end-2);

        dt = 1.0*dx/abs(v);
        [ M, N ] = size(x);
        N = max(M,N);
        y_analytic = u_analitic(x,0);

        u_euler = u_analitic(x,0);
        %initialize lax and some parameters... x
        u_num = u_analitic(x,0);
%         solver = RNFDSolver(N,dx,v/abs(v),abs(v), u_num);
%                 solver = LaxSolver(N,dx,v/abs(v),abs(v), u_num);
        %       solver = HighOrderSolverImplicit(N,dx,dt,direction,velocity, U_0)
        %       solver = HighOrderSolver(N,dx,dt,direction,velocity, U_0)
                solver = SpectralMSolver(N,v/abs(v),a,b,3,abs(v),[],u_num);
        t = dt;
        ctr = 1;
        num_err = 0;
        while(t <= tend)
            %% main loop dx convergence...
            if(mod(ctr,plot_ctr) == 0)
                plot(x,real(u_num),'g',x,real(y_analytic),'--b');
            %    axis([a b -1 1])
            %   legend('lax-wendroff ','euler','analytic sol');
        u_euler = u_analitic(x,0);        
        %initialize lax and some parameters... x
        u_num = u_analitic(x,0);
%         solver = RNFDSolver(N,dx,v/abs(v),abs(v), u_num);
        solver = LaxSolverII(N,dx,v/abs(v),abs(v), u_num);
%       solver = HighOrderSolverImplicit(N,dx,dt,direction,velocity, U_0)
%       solver = HighOrderSolver(N,dx,dt,direction,velocity, U_0)
%         solver = SpectralMSolver(N,v/abs(v),a,b,3,abs(v),u_num);
%         solver = SpectralMSolverGPUKernel(N,v/abs(v),a,b,3,abs(v),u_num);
        t = dt;
        ctr = 1;
        err_euler = 0;
        err_lax = 0;
        compute_time_order = 0.0;
        while(t <= tend)
        %% main loop dx convergence... 
            if(mod(ctr,1000) == 0)
                plot(x,real(u_num),'g',x,real(y_analytic),'--b');
%                 axis([a b -1 1])
%                 legend('lax-wendroff ','euler','analytic sol');

                title(['accuracy comparison: @ t = ' num2str(t)])
                getframe;
            end
            
            num_err(ctr) = 0;
            %analytic solution
            y_analytic = u_analitic(x,t);
%             u_num = solver.make_step(f(x,t-dt),f_t(x,t-dt),0,dt);
                        u_num = solver.make_step(f(x,t-dt), dt);
%                          u_num = solver.set_bdry(y_analytic(1),y_analytic(end));
            if(v < 0)
                u_num = solver.set_bdry('no',y_analytic(end));
            else
                u_num = solver.set_bdry(y_analytic(1),'no');
            end
            %error estimate
            num_err(ctr) =norm(u_num - y_analytic);
=======
            err_lax(ctr) = 0;
            %analytic solution
            y_analytic = u_analitic(x,t);
%             u_num = solver.make_step(f(x,t-dt),f_t(x,t-dt),0,dt);
            tic;
            u_num = solver.make_step(f(x,t-dt), dt);
            %              u_num = solver.set_bdry(y_analytic(1),y_analytic(end));  
%             u_num = solver.set_bdry('no',y_analytic(end));  
            u_num = solver.set_bdry(y_analytic(1),'no');
            compute_time_order = compute_time_order + toc;
            %error estimate
%             err_lax(ctr) = norm(gather((u_num(2:end-1)) - y_analytic(2:end-1)));
>>>>>>> 56594ff310ff69a8566073835808f038b74e84a1
            t = t+dt;
            ctr = ctr+1;
        end
        %store the maximum errors...
<<<<<<< HEAD
        dx_error_samples(order,:) = [max(num_err)];
        
    end
    %     figure;
    %     plot(sizes,dx_error_samples(:,1),'r');
    %     xlabel('N (nr. grid points)' )
    %     ylabel('|\epsilon| - global error');
    %     legend(['\epsilon_{num}'] );
    %     title([' RNFD Convergence with \omega = ' num2str(w)]);
=======
        dx_error_samples(order,:) = [max(err_lax)];
        compute_times(order) = compute_time_order;
    end
%     figure;
%     plot(sizes,dx_error_samples(:,1),'r');
%     xlabel('N (nr. grid points)' )
%     ylabel('|\epsilon| - global error');
%     legend(['\epsilon_{num}'] );
%     title([' RNFD Convergence with \omega = ' num2str(w)]);
>>>>>>> 56594ff310ff69a8566073835808f038b74e84a1
    display(dx_error_samples)
end

%% dt-convergence measurement region
<<<<<<< HEAD
% % % if(t_measure > 0)
% % %     clc;
% % %     %%%% MEASURE THE ERROR IN THE TIME COORDINATE
% % %     orderend = 10;
% % %     dt_init = 10^-1; % initial timestep is 0.009...
% % %     dt_samples(1) = dt_init;
% % %     for i = 2:orderend
% % %         dt_samples(i) = dt_samples(i-1)/2;
% % %     end
% % %     dt_samples = dt_samples';
% % %     
% % %     
% % %     dx = 1.1*dt_init; %%initial dx is 0.0099
% % %     x = 0:dx:xend; x= x'; %transpose x
% % %     dt_error_samples = zeros(orderend,2); %% arrays to store the global error for euler and lax method for each dt configuration
% % %     
% % %     
    %measure the error when we gradually decrease timestep..
% % % %     for order = 1:orderend;
% % % %         %% dt - convergence 1 case
% % % %         clc;
% % % %         display(['time error measurement : order = ' num2str(order)]);
% % % %         dt = dt_samples(order);
% % % %         
% % % %         %boundary condition at left boundary
% % % %         g_l = 0;
% % % %         g_r = @(t) u_analitic(xend,t);
% % % %         %boundary condition at right boundary
% % % %         %initial conditions...
% % % %         u_numeric = zeros(size(x));
% % % %         
% % % %         [ M, N ] = size(x);
% % % %         N = max(M,N);
% % % %         
% % % %         
% % % %         y_analytic = zeros(N,1);
% % % %         u_euler = zeros(N,1);
% % % %         
% % % %         
% % % %         %initialize lax and some parameters... x
% % % %         u_num = zeros(N,1);
% % % %         solver = LaxSolver(N,dx,1,0,v, u_euler);
% % % %         
% % % %         
% % % %         x_short = x(2:N-1);
% % % %         t = 0;
% % % %         ctr = 1;
% % % %         %% end of dt-convergence 1 case setup
% % % %         while(t <= tend)
% % % %             %% dt convergence main loop
% % % %             err_euler(ctr) = 0;
% % % %             err_lax(ctr) = 0;
% % % %             %analytic solution
% % % %             y_analytic = u_analitic(x,t);
% % % %             %numerical solution using euler
% % % %             tmp = f(x,t);
% % % %             u_euler(2:N) = u_euler(2:N) - v*dt/dx*(u_euler(2:N) - u_euler(1:N-1)) + dt*tmp(2:N);
% % % %             u_euler(1) = y_analytic(1);
% % % %             %numerical solution using lax method...
% % % %             df_dt = (f(x_short,t) - f(x_short,t-dt))/dt;
% % % %             df_dx = (f(x_short,t) - f(x(1:N-2),t)/dx);
% % % %             solver.make_step(f(x,t),df_dt,df_dx,dt)
% % % %             u_num = solver.set_bdry(y_analytic(1));
% % % %             if(mod(ctr,plot_ctr) == 0)
% % % %                 plot(x,real(u_num),'--r',x,real(u_euler),'--g',x,real(y_analytic),'b');
% % % %                 legend('lax-wendroff ','euler_sol','analytic sol');
% % % %                 title(['accuracy comparison: @ t = ' num2str(t)])
% % % %                 getframe;
% % % %             end
% % % %             %error estimate
% % % %             err_euler(ctr) = max(abs(u_euler-y_analytic));
% % % %             err_lax(ctr) = max(abs(u_num - y_analytic));
% % % %             t = t+dt;
% % % %             ctr = ctr+1;
% % % %         end
% % % %         %store the maximum errors...
% % % %         dt_error_samples(order,:) = [max(err_euler) max(err_lax)];
% % % %         
% % % %     end
% % % %     
% % % %     figure; % plot in a new figure
% % % %     plot(dt_samples,dt_error_samples(:,1),'r',dt_samples,dt_error_samples(:,2),'g');
% % % %     %     axis([0 dt_samples(end) 0 50]);
% % % %     xlabel('dt (timestep)');
% % % %     ylabel('|\epsilon| - global error');
% % % %     legend(['\epsilon_{euler}'],['\epsilon_{lax}']);
% % % %     title(['Lax-Wendroff scheme convergence dx fixed @ ' num2str(dx) ' ; dt varying']);
% % % %     
% % % % end
=======
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
        solver = LaxSolver(N,dx,1,0,v, u_euler);

        
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
            u_euler(2:N) = u_euler(2:N) - v*dt/dx*(u_euler(2:N) - u_euler(1:N-1)) + dt*tmp(2:N);
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
>>>>>>> 56594ff310ff69a8566073835808f038b74e84a1


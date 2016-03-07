clear; clc;close all;
scenariofile = '..\MB Inputs\standarddeph.set';
simfile = '..\MB Inputs\mlvl_11p0.sim';

amplitudes = [3];
% amplitudes = 5E-5;


amplitudes = [1 ];

taus = {[.4,1.25,1.25]};

compensation = 'constant';
simnames = {[compensation '-disp']};
scen_ctr = 0;
for tau_idx = 1:length(taus)
    for ampl_idx = 1:length(amplitudes)
        
        scen_ctr = scen_ctr+1;
        
        ampl = amplitudes(ampl_idx); tau = taus(tau_idx);
        args = {'init',1};
        lvar = length(args);
        %handle user input!
        if(length(args)>0)
            if mod(lvar,2) ~=0
                display('incorrect number of input arguments. Aborting!');
                return;
            end
            for i = 1:2:lvar
                name = args{i}; val  = args{i+1};
                name = strtrim(lower(name)); %change to lower case!
                if strcmp(name,'init')
                    init = val;
                else if strcmp(name,'workspace')
                        %regexp to load all variables except init!
                        load(val,'-regexp','^(?!.*init.*).*$');
                        
                        
                    else
                        display( ['Unknown input option name: ' name '. Aborting!' ]);
                        return;
                    end
                end
            end
        else
            init = 1;
        end
        
        %parse all input files
        settings = parseSimParams(scenariofile);
        settings = parseSimData(simfile,settings);
        
        
        tch = settings.tch; % characteristic time..
        lch = settings.lch; % characteristic length...
        
        %speed of light in vacuum: 299792458m/s --> in tch/lch
        c_0 = Constants('c',{'time',tch},{'length',lch});
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEVICE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%
        %cavity length  -->
        Ltot = 5; % mm
        
        %phase velocity inside the medium ( in mm per picosecond ... )
        n = settings.n;  c = c_0/n; T_R = 2*Ltot/c; f_R = 1/T_R;
        
        %%%%dipole mtx elements (in Cnm)
        zUL = settings.zUL;
        % hbar in eV-ps
        hbar = Constants('hbar',{'time',tch})/Constants('q0');
        
        %obtain the energies of the core levels for the simulation
        HTB = settings.HTB;
        NLVLS = sqrt(length(HTB)); % nr of levels to consider !
        %reshape into a matrix
        HTB = reshape(HTB,NLVLS,NLVLS).';
        
        if abs(HTB(1,3))>abs(HTB(2,3))
            INJ = 1; IGNORELEVEL = 7; DEPOP = 6;
        else
            INJ =2; IGNORELEVEL = 6; DEPOP = 7;
        end
        
        ULL = 3;
        LLL = 4;
        
        E1 = HTB(INJ,INJ)/hbar;
        E3 = HTB(ULL,ULL)/hbar;
        E2 = HTB(LLL,LLL)/hbar; % in rad/ps
        
        
        %%%% energies: E1 > E3 > E2 Optical transition is 3->2 TUNNELING transition is 1->3
        O13 = HTB(INJ,ULL)/hbar; % in rad/ps
        E13 = E1-E3; %rad/ps; 1->3 traisition freq
        E12 = E1-E2; %rad/ps; 1->2 transition freq
        E32 = E3-E2; %rad/ps  3->2 transition freq (optical transition)
        % central frequency and wave number!
        E0 =  E32; % central OPTICAL frequency.
        
        %gain recovery times
        W = settings.Wmtx;
        NLVLS = sqrt(length(W)); % nr of levels to consider !
        %reshape into a matrix
        W = reshape(W,NLVLS,NLVLS).';
        
        %extract the indices of the rest of the laser levels, which will be
        %included in the system as simple rate equations
        global_idx_rest = setdiff(1:NLVLS,[IGNORELEVEL,INJ,ULL,LLL,DEPOP]);
        N_rest = length(global_idx_rest); %take the number of residual levels!
        
        G = zeros(NLVLS,1);
        % W(DEPOP,:) = zeros(NLVLS,1); W(INJ,DEPOP) = 0;
        G(INJ) = sum(W(INJ,setdiff(1:NLVLS,[INJ,DEPOP,IGNORELEVEL])));
        for lvl = setdiff(1:NLVLS,[INJ,DEPOP,IGNORELEVEL])
            G(lvl) = sum(W(lvl,setdiff(1:NLVLS,[IGNORELEVEL]))); % error here remove the influence of the extra injector level
        end
        
        % % % % Added Pure Dephasing
        if(settings.deph>0)
            gamma_13 = 0.5*(G(INJ)+G(ULL)) +1/tau{scen_ctr}(1); %% dephsing of the resonant tunneling transition
            gamma_32 = 0.5*(G(LLL)+G(ULL))+1/tau{scen_ctr}(2); % dephasing of the optical transision...
            gamma_12 = 0.5*(G(LLL)+G(INJ))+1/tau{scen_ctr}(3); % dephasing of the latest transition
        else
            gamma_13 = 0.5*(G(INJ)+G(ULL)); %% dephsing of the resonant tunneling transition
            gamma_32 = 0.5*(G(LLL)+G(ULL)); % dephasing of the optical transision...
            gamma_12 = 0.5*(G(LLL)+G(INJ)); % dephasing of the latest transition
        end
        
        dE13 = -1i*E13 - gamma_13; %
        dE32 = +1i*(E0 - E32) - gamma_32; %
        dE12 = +1i*(E0 - E12)- gamma_12; %
        
        if (settings.shb > 0)
            shb = 1;
        else
            shb = -1;
        end
        
        Ncarriers = settings.dN*(100^3)*settings.Ld/settings.Lp; % cm^-3 --> m^-3; carrier density
        Overlap = settings.Overlap;  % overlap factor -> dimensionless
        
        gamma = 3*10^-11; % [m/V^2]
        
        alpha = gamma*(Constants('hbar')/( settings.zUL*10^-9*Constants('q0') ))^2;
        alpha = alpha*(10^24/(10^3));%renormalize :)
        
        %calculate the normalization constant, i.e. the number on which we divide
        %all quantities ( except the electric field envelope) to simplify
        %our initial system. it is important as it determines the initial value of
        %the overall electron population inside the system-> a quantity that shall
        %be perserved throughout the whole simulaiton ! !
        
        trace_rho = ((E0*1E12*Ncarriers*Overlap*((zUL*1E-9*Constants('q0'))^2))/(Constants('eps0')*n*Constants('c')*Constants('hbar')))/(1/(lch*tch));
        
        %cavity loss l_0 in (cm^-1) --> l_0*100 in (m^-1) --> 1 mm^-1
        l_0 = settings.loss*100/(1/lch);
        
        %%%%%%%%%%%%%%%%% simuation parameters %%%%%%%%%%%%%%%%%
        
        %grid size in x direction
        N = settings.N; %nr of points in x direction
        x = linspace(0,Ltot,N)';
        dx = x(2) - x(1); dt = dx/c;
        
        k0 = E0/c; D = settings.D*10^2/(1/tch); diffusion = 4*k0^2*D;
        
        %%%% specify some of the mainloop control parameters %%%%
        iter_per_rt = round(T_R/dt);
        checkptIter = iter_per_rt*100; %make a checkpoint every 100RTs.
        tEnd = settings.simRT*T_R; % end time in tps
        plotCtr = settings.plotCtr; %% set on how many iterations should the program plot the intensity
        %simulation info storage arrays -> preallocate
        recordingduration = settings.recordRT*T_R; % how many ps should we record the pulse for
        
        if(init > 0)
            
            
            clc;
            
            x_0 = Ltot/7; 
%             
            tp =0.3;
            aE_in = @(z,time) exp(-(time-(z-x_0)/c).^2/tp^2);
            U = aE_in(x,0)*exp(1i*pi/4);
            V = 0*U;
          
%             U = zeros(N,1);V = U;
%             dw = 2*pi/T_R; dKappa = dw/c;
%             nrmodes = 200; 
%             for mode = -nrmodes/2:nrmodes/2-1
%                 U =  U+ exp(1i*mode*dKappa*(x-x_0));
%             end
            %normalize
            U = ampl*U/max(abs(U));
            
            %             initdata = wgn(N*2-2,1,1);
            %
            %             %sampling rate;
            %             ps = 1/dx;
            %             %number of samples
            %             Ns = length(initdata);
            %             samples = (1:1:Ns)/Ns;
            %             dfilt = fir1(80,[0.001,0.1]);
            %             initdata = filter(dfilt,1,initdata);
            %             U = initdata(1:N);  V = [U(1); initdata(end:-1:N)];
            %
            % population vectors together with first derivative vectors
            r110 = trace_rho*(ones(N,1)*1/3); r110(1) = r110(end);
            r330 = trace_rho*(ones(N,1)*1/3); r330(1) = r330(end);
            r220 = trace_rho*(ones(N,1)*1/3); r220(1)=r220(end);
            
            r11p = zeros(N,1);r33p = zeros(N,1);r22p = zeros(N,1);
            %rate populations
            populations = zeros(N,N_rest);
            gratings = zeros(N,N_rest);
            
            % coherence terms together with first derivative vectors!
            r130 =0*((ones(N,1)-0.5) +1i*(ones(N,1)-0.5)); r130(1) = r130(end);
            r13p = zeros(N,1); r13m = zeros(N,1);
            
            n32p =  1E-15*((rand(N,1)-0.5) +1i*(rand(N,1)-0.5)); n32p(1) = n32p(end);
            n32m =  1E-15*((rand(N,1)-0.5) +1i*(rand(N,1)-0.5)); n32m(1) = n32m(end);
            
            n12p =  1E-15*((rand(N,1)-0.5) +1i*(rand(N,1)-0.5)); n12p(1) = n12p(end);
            n12m =  1E-15*((rand(N,1)-0.5) +1i*(rand(N,1)-0.5)); n12m(1) = n12m(end);
            
            t = dt; nr_steps = settings.nr_steps;
            
            idx = 1; % the index of the predefined point we sample the feild at!
            iter_ctr = 0;    ctr = 1;   % nr of iterations after which to record new statistics!
            %preallocate memory
            
            E_p= 1;  E_m = 1; r11_time = 1;   r33_time = 1;r22_time = 1;
            r11_p_time =1; r33_p_time = 1;  r22_p_time = 1;
            pop_time = cell(N_rest,1);
            for i=1:N_rest
                pop_time{i} =1;
            end
            r130_time = 1;
            r13p_time =1;
            r13m_time = 1;
            
            r110_solver = MS(nr_steps,N,[],r110); r11p_solver =MS(nr_steps,N,[],r11p);
            r330_solver = MS(nr_steps,N,[],r330); r33p_solver = MS(nr_steps,N,[],r33p);
            r220_solver = MS(nr_steps,N,[],r220); r22p_solver = MS(nr_steps,N,[],r22p);
            
            %the solvers for the rest of the auxilliary rates.
            
            rate_eqn_solvers = [];
            grating_eqn_solvers = [];
            for p=1:N_rest
                rate_eqn_solvers{p} =  MS(nr_steps,N,[],populations(:,p));
                grating_eqn_solvers{p} = MS(nr_steps,N,[],gratings(:,p));
            end
            
            n12p_solver = MS(nr_steps,N,[],n12p); n12m_solver = MS(nr_steps,N,[],n12m);
            n32p_solver = MS(nr_steps,N,[],n32p); n32m_solver = MS(nr_steps,N,[],n32m);
            
            r130_solver = MS(nr_steps,N,[],r130);
            r13p_solver = MS(nr_steps,N,[],r13p); r13m_solver = MS(nr_steps,N,[],r13m);
            
            U_solver = RNFDSolver(N,dx,+1,c, U);
            V_solver = RNFDSolver(N,dx,-1,c,V);
            
        end
        
        iterperrecord = 1;
        recordingiter  = round(recordingduration/iterperrecord/dt);
        padsize = recordingiter-length(E_p);
        
        
        %preallocate memory
        E_p= padarray(E_p,padsize,'post');  E_m = padarray(E_m,padsize,'post');
        %store population info
        r11_time = padarray(r11_time,padsize,'post');    r33_time = padarray(r33_time,padsize,'post');
        r22_time = padarray(r22_time,padsize,'post');    r11_p_time = padarray(r11_p_time,padsize,'post');
        r33_p_time = padarray(r33_p_time,padsize,'post');  r22_p_time = padarray(r22_p_time,padsize,'post');
        
        for p = 1:N_rest
            pop_time{p} = padarray(pop_time{p},padsize,'post');
        end
        
        r130_time = padarray(r130_time,padsize,'post');
        r13p_time = padarray(r13p_time,padsize,'post');
        r13m_time = padarray(r13m_time,padsize,'post');
        
        
        transform = @fft;
        itransform = @ifft;
        
        if strcmp(compensation,'perfect')
            compensations = 1;
            compensate = true;
            Yold = transform([U(1:N-1);V(N:-1:2)]); %flip V
            slope = dt*c;
            N2= 2*N-2;
            k_ = 2*pi/N2*[[0:N2/2-1].';[-N2/2:-1].']*1/dx;
            dk = k_(2)-k_(1);
            
          
        else if strcmp(compensation, 'constant');
                N2 = 2*N-2;
                k_ = (2*pi/N2)*[0:N2/2-1,-N2/2:-1].'*1/dx;
                zer_idx = round(N2)/2;

                dk = k_(2)-k_(1); 
                BW = 2/c*2*pi;
                limIdx = k_>0 & k_ < BW;
                kBW = k_(limIdx); 
                
                U_lin_prev = [U(1:N-1);V(N:-1:2)]; %flip V
                Psi_const = 0*k_;
                Ymod = ones(N2,1);
                compensate = false;
                
                dev = 1/c*2*pi; 
                gaussL = exp(-(k_(1:N2/2)).^2/(2*dev^2));
%                 gaussR = flipud(gaussL(1:end));
%                 gausswide = exp(-(k_(1:N2/2)).^2/(5*dev^2));
%                 gausstot = [gaussR;gaussL];
%                 gausstot = [flipud(gausswide);gausswide];
                
            elseif strcmp(compensation,'add-phase')
                beta2 = -0.1 %ps^2/mm;
                alpha2 = -(beta2/2)*c^3;
                N2= 2*N-2;
                k_ = 2*pi/N2*[[0:N2/2-1].';[-N2/2:-1].']*1/dx;
                Psi_2 = dt*alpha2*k_.^2;
                compensate = true;
                
            end
        end
        
        

        info.settings = settings;
        info.cavity = 'FP';
        info.Ltot = Ltot;
        info.N = N;
        info.SIMTYPE = 'WITHOUT DISP. COMP';
        info.l_0 = l_0;
        
        
        recordPhaseIterations = 500;
        while(t< tEnd)
            
            
            if(strcmp(compensation,'constant'))
                if (iter_ctr+1-recordPhaseIterations==0)
                    Y1 = (transform(real(U_lin_prev)));
                    Y2 = (transform(real([U(1:N-1);V(N:-1:2)])));
                    Psi = unwrap(angle(Y2)-angle(Y1));  
                    
                    nconst = recordPhaseIterations; 
                    
                    psi0 = Psi(1); slope = (Psi(2)-Psi(1))/dk;
                    Psi_const = (Psi(1:N2/2)-slope*k_(1:N2/2)-psi0).*gaussL; 
                    
                    Ymod = fftshift(exp(-1i*[-flipud(Psi_const);Psi_const]/nconst));  % Specify the filter response at the freqs in f1.

%                     F1= linspace(-1,1,N2).'; % Define the normalized frequencies' array
%                     f = fdesign.arbmagnphase('N,F,H',100,F1,H1);
%                     Hd = design(f,'equiripple');
%                     hfvt = fvtool(Hd, 'Color','w');
%                     
%                     hfvt(2) = fvtool(Hd,'Analysis','phase','Color','white');
% 
%                     ax = hfvt(2).CurrentAxes;
%                     ax.NextPlot = 'add';
%                     pidx = find(F1>=0);
%                     plot(ax,F1,[flipud(unwrap(angle(H1(pidx-1:-1:1,1))));unwrap(angle(H1(pidx:end,1)))],'k--')
           
%                     U_lin_prev = [U(1:N-1);V(N:-1:2)]; %flip V
                    compensate = true; 
                end
            end
            
            
            
            %%plot some of the results if neeed ariseth :D
            if(mod(iter_ctr, 1000) == 0)
                clc;
                info.iter_ctr = iter_ctr;
                info.RT = t/T_R;
                intensity = U.*conj(U) + V.*conj(V) ;
                info.maxInt  =  max(intensity);
                printINFO(info);
                
                %                 plotyy(x,[real(U),real(V)],x,[r110,r330,r220,populations]);
                subplot(2,1,1);
                plot(real([U(1:N-1);V(N:-1:2)]))
                subplot(2,1,2);
                plot(x,[(imag(U)),imag(V)]);
                title([ settings.name ' (MS + RNFD)  @ t = ' num2str(t)]);
                getframe;
            end
            
            %%%% obtain the field, field intensity and the total population at position "idx" ...
            
            if((t >= tEnd - recordingduration) && mod(iter_ctr,iterperrecord) == 0)
                
                %store fields info
                E_p(ctr)= U(idx);  E_m(ctr)= V(idx);
                %store population info
                r11_time(ctr)= r110(idx);
                r33_time(ctr)= r330(idx);
                r22_time(ctr) = r220(idx);
                
                r11_p_time(ctr) = r11p(idx);
                r33_p_time(ctr) = r33p(idx);
                r22_p_time(ctr) = r22p(idx);
                
                %calculate total populations
                for p = 1:N_rest
                    pop_time{p}(ctr) = populations(idx,p);
                end
                
                r130_time(ctr) = r130(idx);
                r13p_time(ctr) = r13p(idx);
                r13m_time(ctr) = r13m(idx);
                
                ctr = ctr+1;
            end
            
            %%%%%% BLOCH PART %%%%%%
            %%%% POPULATIONS
            r110_t = 1i*O13.*(r130-conj(r130)) +( W(ULL,INJ) + W(ULL,DEPOP) )*r330 + ( W(LLL,INJ)+W(LLL,DEPOP) )*r220 - G(INJ)*r110;
            for p = 1:N_rest
                p_glob_idx = global_idx_rest(p);
                r110_t = r110_t + (W(p_glob_idx,INJ)+W(p_glob_idx,DEPOP))*populations(:,p);
            end
            r110_solver.make_step(r110_t,dt);
            
            lmInteraction = conj(U).*n32p + conj(V).*n32m;
            r330_t = 1i*O13.*(conj(r130) - r130) +1i/2.*(lmInteraction-conj(lmInteraction)) + r110*W(INJ,ULL) + r220*W(LLL,ULL) - G(ULL)*r330;
            for p = 1:N_rest
                p_glob_idx = global_idx_rest(p);
                r330_t = r330_t + W(p_glob_idx,ULL)*populations(:,p);
            end
            r330_solver.make_step(r330_t,dt);
            
            r220_t = -1i/2.*(lmInteraction-conj(lmInteraction)) + r110*W(INJ,LLL) + r330*W(ULL,LLL) - G(LLL)*r220;
            for p = 1:N_rest
                p_glob_idx = global_idx_rest(p);
                r220_t = r220_t + W(p_glob_idx,LLL)*populations(:,p);
            end
            r220_solver.make_step(r220_t,dt);
            
            for p = 1:N_rest
                p_glob_idx =global_idx_rest(p);
                rpp_0_t = W(INJ,p_glob_idx)*r110+W(ULL,p_glob_idx)*r330+W(LLL,p_glob_idx)*r220 - G(p_glob_idx)*populations(:,p);
                %add the rate equations part!
                for j =1:N_rest % nr
                    if j~= p
                        j_glob_idx =global_idx_rest(j);
                        rpp_0_t = rpp_0_t + W(j_glob_idx,p_glob_idx)*populations(:,j);
                    end
                end
                rate_eqn_solvers{p}.make_step(rpp_0_t,dt);
            end
            
            %%%% COHERENCES
            r130_t = dE13*r130 + 1i*O13*(r110-r330) +1i/2*(conj(U).*n12p + conj(V).*n12m);
            r130_solver.make_step(r130_t,dt);
            
            n32p_t = dE32*n32p + 1i/2*(U.*(r330-r220) + V.*(r33p-r22p)) - 1i*O13*n12p;
            n32p_solver.make_step(n32p_t,dt);
            
            n32m_t = dE32*n32m + 1i/2*(V.*(r330-r220) + U.*conj(r33p-r22p)) - 1i*O13*n12m;
            n32m_solver.make_step(n32m_t,dt);
            
            n12p_t = dE12*n12p +1i/2*(U.*r130 + V.*r13p) - 1i*O13*n32p;
            n12p_solver.make_step(n12p_t,dt);
            
            n12m_t = dE12*n12m + 1i/2*(V.*r130 +U.*r13m) - 1i*O13*n32m;
            n12m_solver.make_step(n12m_t,dt);
            
            if(settings.shb > 0 )
                
                %%% r11+
                r11p_t = 1i*O13.*(r13p-conj(r13m)) + (W(ULL,INJ)+W(ULL,DEPOP))*r33p+ (W(LLL,INJ)+W(LLL,DEPOP))*r22p - (G(INJ)+diffusion)*r11p;
                r11p_solver.make_step(r11p_t,dt);
                
                %%% r33+
                r33p_t = 1i*O13.*(conj(r13m)-r13p)+1i/2*(conj(V).*(n32p) -(U).*conj(n32m)) +  W(INJ,ULL)*r11p + W(LLL,ULL)*r22p - (G(ULL)+diffusion)*r33p;
                r33p_solver.make_step(r33p_t,dt);
                
                %%% r22+
                r22p_t = -1i/2*(conj(V).*(n32p) -(U).*conj(n32m)) +  W(INJ,LLL)*r11p + W(ULL,LLL)*r33p - (G(LLL)+diffusion)*r22p;
                r22p_solver.make_step(r22p_t,dt);
                
                %%% r13+
                r13p_t = (dE13-diffusion)*r13p + 1i*O13*(r11p-r33p) +1i/2*(conj(V).*n12p);
                r13p_solver.make_step(r13p_t,dt);
                
                %%% r13-
                r13m_t = (dE13-diffusion)*r13m + 1i*O13*conj(r11p-r33p) + 1i/2*(conj(U).*n12m);
                r13m_solver.make_step(r13m_t,dt);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            %              I =  abs(U).^2+abs(V).^2;
            %             l = (l_0  - alpha*I);
            %             l = l_0/(1+I/I_sat);
            losses = -c*l_0*ones(N,1);
            %             losses = -c*l;
            
            
            U = U_solver.make_step(-1i*c*n32p,-1i*c*n32p_t,losses,dt);
            V = V_solver.make_step(-1i*c*n32m,-1i*c*n32m_t,losses,dt);
            
            %set the boundaries... and obtain the final solution
            U = U_solver.set_bdry(V(1),'no');
            V = V_solver.set_bdry('no',U(N));
            %%%%%%%%%%%%%%%%%%%%%%%%
            
            if (compensate)
                if(strcmp(compensation,'perfect')) %perfect phase compensation
                    Ynew = transform([U(1:N-1);V(N:-1:2)]);
                    Ynew = (abs(Ynew).*exp(1i*slope*k_+1i*unwrap(angle(Yold))));
                    U_lin = itransform(Ynew);
                    U = U_lin(1:N);  V = [U(1); U_lin(end:-1:N)];
                    U_solver.set_latest_solution(U);
                    V_solver.set_latest_solution(V);
                    Yold = Ynew;
                    
                elseif strcmp(compensation,'add-phase') % add quadratic phase
                    Ynew = transform([U(1:N-1);V(N:-1:2)]);
                    Ynew = Ynew.*exp(1i*Psi_2);
                    U_lin = itransform(Ynew);
                    U = U_lin(1:N);  V = [U(1); U_lin(end:-1:N)];
                    U_solver.set_latest_solution(U);
                    V_solver.set_latest_solution(V);
                else %constant phase comp
                    U_lin = itransform(transform([U(1:N-1);V(N:-1:2)]).*Ymod);
%                     U_lin = filter(Hd,[U(1:N-1);V(N:-1:2)]);
                    U = U_lin(1:N);  V = [U(1); U_lin(end:-1:N)];
                    U_solver.set_latest_solution(U);
                    V_solver.set_latest_solution(V);

                end
            end
            
            
            r110 = r110_solver.get_latest_solution();
            r330 = r330_solver.get_latest_solution();
            r220 = r220_solver.get_latest_solution();
            
            for p = 1:length(global_idx_rest)
                populations(:,p) = rate_eqn_solvers{p}.get_latest_solution();
            end
            
            r130 = r130_solver.get_latest_solution(); n32p = n32p_solver.get_latest_solution();
            n32m = n32m_solver.get_latest_solution(); n12p = n12p_solver.get_latest_solution(); n12m = n12m_solver.get_latest_solution();
            
            if(settings.shb > 0)
                r11p = r11p_solver.get_latest_solution();
                r33p = r33p_solver.get_latest_solution();
                r22p = r22p_solver.get_latest_solution();
                r13p = r13p_solver.get_latest_solution();
                r13m = r13m_solver.get_latest_solution();
            end
            
            t = t+dt;
            iter_ctr = iter_ctr + 1;
            
        end
        savename = [settings.name '_' simnames '_N_' num2str(settings.N) '_FP_' num2str(settings.simRT) ];
        save(savename);
            
    end
end


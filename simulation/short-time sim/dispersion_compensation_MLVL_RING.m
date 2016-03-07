clear; clc;close all;
scenariofile = '..\MB Inputs\standarddeph.set';
simfile = '..\MB Inputs\mlvl_11p0.sim';

amplitudes = [ 0.005 0.05 0.5 1 2 3 4 5 6 ];
% amplitudes = 5E-5;

taus = {[.4,1.25,1.25]};

compensation = 'add-phase';
simnames = {[compensation '-compensation']};
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
                name = args{i};
                val  = args{i+1};
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
        Ltot = 10; % mm
        
        
        %phase velocity inside the medium ( in mm per picosecond ... )
        n = settings.n; c = c_0/n; T_R = Ltot/c; f_R = 1/T_R;
        
        %%%%dipole mtx elements (in Cnm)
        zUL = settings.zUL;
        % hbar in eV-ps
        hbar = Constants('hbar',{'time',tch})/Constants('q0');
        
        HTB = settings.HTB;
        NLVLS = sqrt(length(HTB)); % nr of levels to consider !
        %reshape into a matrix
        HTB = reshape(HTB,NLVLS,NLVLS).';
        
        if abs(HTB(1,3))>abs(HTB(2,3))
            INJ = 1; DEPOP = 6; IGNORELEVEL = 7;
        else
            INJ =2; DEPOP = 7; IGNORELEVEL = 6;
        end
        ULL = 3; LLL = 4;
        
        E1 = HTB(INJ,INJ)/hbar;
        E3 = HTB(ULL,ULL)/hbar;
        E2 = HTB(LLL,LLL)/hbar; % in rad/ps
        
        %%%% energies: E1 > E3 > E2 Optical transition is 3->2 TUNNELING transition is 1->3
        O13 = HTB(INJ,ULL)/hbar; % in rad/ps
        E13 = E1-E3; %rad/ps; 1->3 traisition freq
        E12 = E1-E2; %rad/ps; 1->2 transition freq
        E32 = E3-E2; %rad/ps  3->2 transition freq (optical transition)
        % central frequency and wave number!
        E0 = E32; % central OPTICAL frequency.
        
        %gain recovery times
        W = settings.Wmtx;
        NLVLS = sqrt(length(W)); % nr of levels to consider !
        %reshape into a matrix
        W = reshape(W,NLVLS,NLVLS).';
        
        %extract the indices of the rest of the laser levels, which will be
        %included in the system as simple rate equations
        idx_rest = setdiff(1:NLVLS,[INJ,ULL,LLL,DEPOP,IGNORELEVEL]);
        N_rest = length(idx_rest); %take the number of residual levels!
        
        G = zeros(NLVLS,1);
        G(INJ) = sum(W(INJ,setdiff(1:NLVLS,[INJ,DEPOP,IGNORELEVEL])));
        for lvl = setdiff(1:NLVLS,[INJ,DEPOP,IGNORELEVEL])
            G(lvl) = sum(W(lvl,setdiff(1:NLVLS,[IGNORELEVEL])));
        end
        
        
        Tdeph_1 = taus{tau_idx}(1);
        Tdeph_2 = taus{tau_idx}(2);
        Tdeph_3 = taus{tau_idx}(3);
        
        % % % % Added Pure Dephasing
        gamma_13 = 1/Tdeph_1; %% dephsing of the resonant tunneling transition
        gamma_32 = 1/2*(G(ULL) + G(LLL)) + 1/Tdeph_2; % dephasing of the optical transision...
        gamma_12 = 1/2*(G(INJ) + G(LLL)) + 1/Tdeph_3; % dephasing of the latest transition
        
        dE13 = -1i*E13 - gamma_13; %
        dE32 = +1i*(E0 - E32) - gamma_32; %
        dE12 = +1i*(E0 - E12)- gamma_12; %
        
        % in the new system of units (mm/picosecond)
        
        Ncarriers = settings.dN*(100^3)*settings.Ld/settings.Lp; % cm^-3 --> m^-3; carrier density
        Overlap = settings.Overlap;  % overlap factor -> dimensionless
        
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
        N = 5000; %nr of points in x direction
        x = linspace(0,Ltot,N)';
        dx = x(2) - x(1); dt = dx/c;
        
        k0 = E0/c;
        
        %%%% specify some of the mainloop control parameters %%%%
        iter_per_rt = round(T_R/dt);
        
        NRT = 5;
        tEnd = NRT*T_R; % end time in tps
        plotCtr = settings.plotCtr; %% set on how many iterations should the program plot the intensity
        pulse_len = tEnd; % how many ps should we record the pulse for
        
        if(init > 0)
            
            x_0 = Ltot/7; tp =.6;
            aE_in = @(z,time) exp(-(time-(z-x_0)/c).^2/tp^2);
            U = aE_in(x,0);
            U = ampl*U/max(abs(U));
            U_uncomp = U;
            U_prev = U;
            
            U_t = zeros(5,1);  V_t = zeros(5,1);
            checkpoints = [1, round(N/5), round(2*N/5),  round(2.5*N/5),round(3.5*N/5)];
            
            Yx = 0*x; Yx(checkpoints) = 1;
            
            ctr = 0;    single_ctr = 1; idx = 1;
            
            % population vectors together with first derivative vectors
            r33 = trace_rho*0.4*ones(N,1); r33(1) = r33(end);
            r22 = trace_rho*0.2*ones(N,1); r22(1) = r22(end);
            r11 = trace_rho*0.4*ones(N,1); r11(1)=r11(end);
            
            %rate populations
            populations = zeros(N,N_rest);
            
            % coherence terms
            r13 = 1E-20*((rand(N,1)-0.5) +1i*(rand(N,1)-0.5)); r13(1) = r13(end);
            n32 = 1E-20*((rand(N,1)-0.5) +1i*(rand(N,1)-0.5)); n32(1) = n32(end);
            n12 = 1E-20*((ones(N,1)-0.5) +1i*(ones(N,1)-0.5)); n12(1) = n12(end);
            
            t = 0; nr_steps = 5;
            
            iter_ctr = 0;
            iterperrecord = 1 ; % nr of iterations after which to record new statistics!
            
            r11_solver = MS(nr_steps,N,[],r11); r22_solver = MS(nr_steps,N,[],r22);
            r33_solver = MS(nr_steps,N,[],r33); n12_solver = MS(nr_steps,N,[],n12);
            r13_solver = MS(nr_steps,N,[],r13); n32_solver = MS(nr_steps,N,[],n32);
            
            rate_eqn_solvers = [];
            for p=1:N_rest
                rate_eqn_solvers{p} =  MS(nr_steps,N,[],populations(:,p));
            end
            
            wave_solver = RNFDSolver(N,dx,c/abs(c),abs(c), U);
            uncompensated_wave_solver = RNFDSolver(N,dx,c/abs(c),abs(c), U);
            
            
        end
        
        
        probe_pt  = 5 ;
        
        transform = @fft;
        itransform = @ifft;
        
        if strcmp(compensation,'perfect')
            compensations = 1;
            slope = -dt*c;
            compensate = true;
            oldphase = unwrap(angle(transform(U)));
            k_ = 2*pi/N*[[0:N/2-1].';[-N/2:-1].']*1/dx;
            dk = k_(2)-k_(1);
            
        else if strcmp(compensation, 'constant');
                recordPhaseIterations = 1;
                load(['Phase_const']);
                compensate = true;
            elseif strcmp(compensation,'add-phase')
                
                %init k - postivie and then negative freq. 
                k_ = 2*pi/N*[[0:N/2-1].';[-N/2:-1].']*1/dx;
               
                
                %setup the windowing func
                boxcaridx = fftshift(k_*c/2/pi<50 & k_*c/2/pi >-50);
                boxbox = boxcaridx(boxcaridx);
                win = 0*k_;
                win(boxcaridx) = fftshift(boxbox.*hanning(length(boxbox)));
                
                %calculate the disperson coeff and set the dispersion phase
                beta2 = -1;%ps^2/mm;
                alpha2 = -(beta2/2)*c^3;
                %                 Psi_L = (k_(1:N/2).^2).*alpha2;
                %                 Psi_2 =[Psi_L;-flipud(Psi_L)]*dt;
                Psi_2 = (k_.^2).*alpha2*dt.*win;
                
                Ymod= win.*exp(1i*Psi_2);
                compensate = true;
                
                
            else
                k_ = (2*pi/N)*[0:N/2-1,-N/2:-1].'*1/dx;
                zer_idx = round(N)/2;
                dk = k_(2)-k_(1);
                
                boxcaridx = fftshift(k_*c/2/pi<50 & k_*c/2/pi >-50);
                boxbox = boxcaridx(boxcaridx);
                win = 0*k_;
                win(boxcaridx) = fftshift(boxbox.*hanning(length(boxbox)));
                
                Psi_const = 0*k_;
                Ymod = ones(N,1);
                compensate = false;
                dev = 2/c*2*pi;
                
            end
        end
        
        info.settings = settings;
        info.cavity = 'RING';
        info.Ltot = Ltot;
        info.N = N;
        info.l_0 = l_0;
        info.SIMTYPE = ['With ' compensation ' phase compensation. '];
        
        
        recordPhaseIterations = 6000;
        startrecordPhase = 1;
        
        while(t< tEnd)
          
            
            if(~compensate)
                  
                if (iter_ctr == startrecordPhase)
                    U_prev = U;
                end

                if( iter_ctr-startrecordPhase == recordPhaseIterations )
                    % % %                 old version -tested and supposedly working starts here
                    % % %                 phase compensation is designed for a real signal
                    % % %                 i.e. it posseses the symmetry properties
                    % that a real linear finite impulse response filter ought to have. It is non-causeal tho.
                    %
                    %                     Y1 = transform(real(U_prev)); Y2 = transform(real(U));
                    %                     Psi = unwrap(angle(Y2)-angle(Y1));
                    %                     nconst = recordPhaseIterations;
                    %                     psi0 = Psi(1); slope = (Psi(2)-Psi(1))/dk;
                    %                     Psi_const_ug = (Psi(1:N/2)-slope*k_(1:N/2)-psi0);
                    %                     Psi_const = (Psi_const_ug).*gaussL;
                    %
                    %                     Psi = [-flipud(Psi_const);Psi_const]/nconst;
                    %                     Psi_ug = [-flipud(Psi_const_ug);Psi_const_ug]/nconst;
                    %
                    %                     Ymod = fftshift(exp(-1i*Psi));  % Specify the filter response at the freqs in fconst','Ymod','k_');
                    %                     Ymod = fftshift(win.*exp(-1i*Psi_ug));
                    %
                    %
                    %                     F1=linspace(-1,1,N).'; % Define the first stopband frequencies.
                    %                     H1 = [gaussR;gaussL].*exp(-1i*Psi*nconst);  % Specify the filter response at the freqs in f1.
                    %                     f = fdesign.arbmagnphase('N,F,H',120,F1,H1);
                    %                     Hd = design(f,'equiripple');
                    %                     hfvt = fvtool(Hd, 'Color','w');
                    %
                    %                     hfvt(2) = fvtool(Hd,'Analysis','phase','Color','white');
                    %
                    %                     ax = hfvt(2).CurrentAxes;
                    %                     ax.NextPlot = 'add';
                    %                     pidx = find(F1>=0);
                    %                     plot(ax,F1,[flipud(unwrap(angle(H1(pidx-1:-1:1))));unwrap(angle(H1(pidx:end)))],'k--')
                    %                     save(['Phase_const' strrep(num2str(ampl),'.','p') '.mat'],'Ymod','Psi','Psi_ug','Hd');
                    
                    
                    Y1 = transform(real(U_prev)); Y2 = transform(real(U));
                    Psi = unwrap(fftshift(angle(Y2)-angle(Y1)));
                    Psi = fftshift(Psi);
                    psi0 = Psi(1); slope = (Psi(2)-Psi(1))/dk;
                    Psi_ = Psi - psi0 - slope*k_;
                    Ymod = exp(-1i*Psi_.*win/recordPhaseIterations).*win;
                    save(['Phase_const' strrep(num2str(ampl),'.','p') '.mat'],'Ymod','Psi_');
                    t = tEnd;
%                     return;
                end
            end
            
            %%plot some of the results if neeed ariseth :D
            if(mod(iter_ctr,100) == 0)
                clc;
                info.iter_ctr = iter_ctr;
                info.RT = t/T_R;
                intensity = U.*conj(U) ;
                info.maxInt  =  max(intensity);
                printINFO(info);
                
                plotyy(x,[real(U),imag(U),real(U_uncomp)],x,Yx);
                title([ settings.name ' (MS + RNFD)  @ t = ' num2str(t)]);
                getframe;
            end
            
            if(t >= tEnd - pulse_len && mod(iter_ctr,iterperrecord) == 0)
                
                E_p(single_ctr) = U(checkpoints(probe_pt));
                %store population info
                r11_time(single_ctr)= r11(checkpoints(probe_pt));
                r33_time(single_ctr)= r33(checkpoints(probe_pt));
                r22_time(single_ctr) = r22(checkpoints(probe_pt));
                
                %calculate total populations
                for p = 1:N_rest
                    pop_time(p,single_ctr) = populations(checkpoints(probe_pt),p);
                end
                
                single_ctr = single_ctr+1;
            end
            
            
            
            U_t(1,ctr+1) = U(checkpoints(1));
            U_t(2,ctr+1) = U(checkpoints(2));
            U_t(3,ctr+1) = U(checkpoints(3));
            U_t(4,ctr+1) = U(checkpoints(4));
            U_t(5,ctr+1) = U(checkpoints(5));
            ctr = ctr  +1 ;
            
            %%%%%% BLOCH PART %%%%%%
            %%%% POPULATIONS
            r11_t = 1i*O13*(r13-conj(r13)) +( W(ULL,INJ) + W(ULL,DEPOP) )*r33 + ( W(LLL,INJ)+W(LLL,DEPOP) )*r22 - G(INJ)*r11;
            for p = 1:N_rest
                p_glob_idx = idx_rest(p);
                r11_t = r11_t + (W(p_glob_idx,INJ)+W(p_glob_idx,DEPOP))*populations(:,p);
            end
            r11_solver.make_step(r11_t,dt);
            
            lmInteraction = conj(U).*n32;
            r33_t = -1i*O13*(r13-conj(r13)) +1i/2.*(lmInteraction-conj(lmInteraction)) + r11*W(INJ,ULL) + r22*W(LLL,ULL) - G(ULL)*r33;
            for p = 1:N_rest
                p_glob_idx = idx_rest(p);
                r33_t = r33_t + W(p_glob_idx,ULL)*populations(:,p);
            end
            r33_solver.make_step(r33_t,dt);
            
            r22_t = -1i/2.*(lmInteraction-conj(lmInteraction)) + r11*W(INJ,LLL) + r33*W(ULL,LLL) - G(LLL)*r22;
            for p = 1:N_rest
                p_glob_idx = idx_rest(p);
                r22_t = r22_t + W(p_glob_idx,LLL)*populations(:,p);
            end
            r22_solver.make_step(r22_t,dt);
            
            for p = 1:N_rest
                
                p_glob_idx =idx_rest(p);
                rpp_0_t = W(INJ,p_glob_idx)*r11+W(ULL,p_glob_idx)*r33+W(LLL,p_glob_idx)*r22 - G(p_glob_idx)*populations(:,p);
                
                %add the rate equations part!
                for j =1:N_rest % nr
                    if j~= p
                        j_glob_idx =idx_rest(j);
                        rpp_0_t = rpp_0_t + W(j_glob_idx,p_glob_idx)*populations(:,j);
                    end
                end
                rate_eqn_solvers{p}.make_step(rpp_0_t,dt);
                
            end
            
            %%% coherences! r13 n32 n12
            r13_t = dE13*r13 + 1i*O13*(r11 - r33) + 1i/2*conj(U).*n12;
            r13_solver.make_step(r13_t,dt);
            
            n32_t = dE32*n32 + 1i/2*U.*(r33-r22) - 1i*O13*n12;
            n32_solver.make_step(n32_t,dt);
            
            n12_t = dE12*n12 + 1i/2*U.*r13 - 1i*O13*n32;
            n12_solver.make_step(n12_t,dt);
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            
            losses = -c*l_0*ones(N,1);
            
            
            %nondispersive wave!
            U_uncomp = uncompensated_wave_solver.make_step(-1i*c*n32*0,-1i*c*n32_t*0,losses*0,dt);
            U_uncomp = uncompensated_wave_solver.set_bdry(U_uncomp(end),'no');
            
            %dispersive wave!
            U = wave_solver.make_step(-1i*c*n32,-1i*c*n32_t,losses,dt);
            U = wave_solver.set_bdry(U(end),'no');
            
            % constant and perferct dispersion compensation
            if (compensate)
                if(strcmp(compensation,'perfect')) %linear compensationa nd the old phase
                    Ynew  = transform(real(U));
                    Ynew = abs(Ynew).*exp(-1i*c*t*k_+1i*(0*oldphase));
                    U = itransform(Ynew);
                    oldphase = unwrap(angle(Ynew));
                elseif strcmp(compensation,'add-phase')
                    %add a constant phase
                    U = itransform(transform(U).*Ymod);
                    
                else
                    U = (itransform(transform(U).*Ymod));
                end
                wave_solver.set_latest_solution(U);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            
            r11 = r11_solver.get_latest_solution();
            r33 = r33_solver.get_latest_solution();
            r22 = r22_solver.get_latest_solution();
            
            for p = 1:length(idx_rest)
                populations(:,p) = rate_eqn_solvers{p}.get_latest_solution();
            end
            
            r13 = r13_solver.get_latest_solution();
            n32 = n32_solver.get_latest_solution();
            n12 = n12_solver.get_latest_solution();
            
            t = t+dt;
            iter_ctr = iter_ctr + 1;
        end
        
        %         mkdir(simnames{tau_idx});
        %         savename = [settings.name '_DA_' simnames{tau_idx} '_' strrep(num2str(ampl),'.','p') '_RING_' num2str(N) '_' strrep(num2str(Tdeph_1),'.','p') '_'  strrep(num2str(Tdeph_2),'.','p') '_'  strrep(num2str(Tdeph_3),'.','p')];
        %         save([simnames{tau_idx} '\' savename]);
        
    end
end
clear; clc;close all;
scenariofile = 'MB Inputs\verystrongdeph_longcavity.set';
simfile = 'MB Inputs\mlvl_11p0.sim';

amplitudes = [.2];
% taus = {[0.175,1.25,1.25],[100,1.25,1.25],[0.4,0.25,1.25],[0.4,100,1.25],[0.4,1.25,0.25],[0.4,1.25,100]};
% simnames = {'1-3-incoherent','1-3-max-coherent','3-2-incoherent','3-2-max-coherent','1-2-incoherent','1-2-max-coherent'};

taus = {[0.4,0.8,0.8],};
simnames = {'disp-comp'};
scen_ctr = 0;
plotON = false;

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
        gamma_13 = 1/2*(G(INJ) + G(ULL)) + 1/Tdeph_1; %% dephsing of the resonant tunneling transition
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
        N = 4000; %nr of points in x direction
        x = linspace(0,Ltot,N)';
        dx = x(2) - x(1); dt = dx/c;
        
        k0 = E0/c;
        
        %%%% specify some of the mainloop control parameters %%%%
        iter_per_rt = round(T_R/dt);
        
        NRT = 10;
        tEnd = NRT*T_R; % end time in tps
        plotCtr = 1000; %% set on how many iterations should the program plot the intensity
        pulse_len = tEnd; % how many ps should we record the pulse for
    
        if(init > 0)
            
            x_0 = Ltot/7; tp =2;
            aE_in = @(z,time) exp(-(time-(z-x_0)/c).^2/tp^2);
            U = zeros(N,1);
            dw = 2*pi/T_R; sig = 2;
            %             nmods = 200;
            %             for m = -nmods/2:nmods/2-1
            %                 U = U+exp(1i*m*dw/c*(x-x_0));
            %             end
            U = aE_in(x,0);
            %normalize
            U = ampl*U/max(abs(U));
            U_uncomp = U;
            U_prev = U;
            k_ = 2*pi*linspace(-1/2,1/2,length(U))'*1/dx;
            dk = k_(2)-k_(1);
            
            
            U_t = zeros(5,1);  V_t = zeros(5,1);
            checkpoints = [1, round(N/5)+1, 2*round(N/5) , 3*round(N/5), 4*round(N/5)];
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
        
        legendinfo = {'real(U)','real(U_{uncomp})','inj','ull','lll'}; nlegend = length(legendinfo);
        for p = 1:length(idx_rest)
            legendinfo{nlegend+p} = [' lvl: ' num2str(idx_rest(p)) ];
        end
        
        probe_pt  = 5 ;
        
        BW = 1/c; %1/mm
        lim_idx= and(k_/2/pi>=-BW/2,k_/2/pi<=BW/2);
        zer_idx = and(k_>=(-dk),k_<=dk);
        U_prev = U;
        Psi_const = 0*k_; 
        
        compensations = 1;
        recordPhaseIterations = 1000; 
        
        info.settings = settings;
        info.cavity = 'RING';
        info.Ltot = Ltot;
        info.N = N;
        info.l_0 = l_0;
        info.SIMTYPE = ['WITH PHASE COMPENSATION every ' num2str(compensations) 'iters'];
        slope = compensations*dt*c;
        transform = @ifft; 
        itransform = @fft;
        
       
    
        while(t< tEnd)
            
            
            if(iter_ctr == recordPhaseIterations)
                    
                Y1 = fftshift(transform(U_prev));
                Y2 = fftshift(transform(U));
                Psi_const = unwrap(angle(Y2))-unwrap(angle(Y1));
                psi0 = mean(Psi_const(zer_idx)); slope = diff(Psi_const(lim_idx))/dk;
                slope = slope(round(length(slope)/2));
                Psi_const = fftshift((Psi_const-slope*k_-psi0).*hanning(length(k_))/recordPhaseIterations); 
            end
          
            
            
%             % dispersion compensation
%             if (mod(iter_ctr+1+recordPhaseIterations,compensations) == 0)
%                 
%                 Y1 = fftshift(transform(U_prev));
%                 Y2 = fftshift(transform(U));
%                 Psi_k = unwrap(angle(Y2))-unwrap(angle(Y1));
%                 psi0 = mean(Psi_k(zer_idx));
%                 
% %               slopes = [(diff(Psi_k)/dk); 0];
%                 %                 slope = slopes(round(length(slopes)/2));
%                 %                 Psi_k(1:end-1)= Psi_k(1:end-1) - slopes.*k_(1:end-1)-psi0;
%                 %                 Psi_k(end) = Psi_k(1);
% %                 slope = slopes(round(length(slopes)/2));
%                 %                 Psi_k(lim_idx)= Psi_k(lim_idx) -[slopes(lim_idx)].*k_(lim_idx);
%                 Psi_k= Psi_k -slope.*k_-psi0;
%                 
%                 %                 Psi_k = Psi_k ;
%                 %                 Psi_k(end) = Psi_k(1);
%                 
%                 Psi_k_comp = Psi_k.*hanning(length(Psi_k));
%                 %                 Psi_k_comp = Psi_k;
% %                 U_prev = itransform(fftshift((Y2.*exp(-1i*Psi_k_comp))));             
%                 U_prev = itransform(fftshift(abs(Y2).*exp(1i*slope*k_+1i*unwrap(angle(Y1)))));
% 
%                 U = U_prev;
%                 wave_solver.set_latest_solution(U);
%             end
            
            %%plot some of the results if neeed ariseth :D
            if(mod(iter_ctr,100) == 0)
                clc;
                clc;
                info.iter_ctr = iter_ctr;
                info.RT = t/T_R;
                intensity = U.*conj(U) ;
                info.maxInt  =  max(intensity);
                printINFO(info);
                
                plotyy(x,[real(U),real(U_uncomp)],x,[r11,r33,r22,populations]);
                legend(legendinfo);
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
            
            
            
            %             rt = floor(t/T_R)+1;
            rt = 1;
            U_t(1,ctr(rt)+1) = U(checkpoints(1));
            U_t(2,ctr(rt)+1) = U(checkpoints(2));
            U_t(3,ctr(rt)+1) = U(checkpoints(3));
            U_t(4,ctr(rt)+1) = U(checkpoints(4));
            U_t(5,ctr(rt)+1) = U(checkpoints(5));
            
            %Update the round-trip counter!
            ctr(rt) = ctr(rt)+1;
            
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
            %
            %dispersive wave!
            U = wave_solver.make_step(-1i*c*n32,-1i*c*n32_t,losses,dt);
            U = wave_solver.set_bdry(U(end),'no');
            
            % constant dispersion compensation
            if (iter_ctr>recordPhaseIterations )               
                U = itransform(transform(U).*exp(-1i*Psi_const));
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
        
        %         dispanalysis_IV;
        
        mkdir(simnames{tau_idx});
        savename = [settings.name '_DA_' simnames{tau_idx} '_' strrep(num2str(ampl),'.','p') '_RING_' num2str(N) '_' strrep(num2str(Tdeph_1),'.','p') '_'  strrep(num2str(Tdeph_2),'.','p') '_'  strrep(num2str(Tdeph_3),'.','p')];
        save([simnames{tau_idx} '\' savename]);
        
    end
end
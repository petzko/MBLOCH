clear; clc;close all;
scenariofile = '..\MB Inputs\verystrongdeph_longcavity.set';
simfile = '..\MB Inputs\mlvl_11p0.sim';

amplitudes = [ 0.005 0.05 0.5 1 2 3 ];
taus = {[0.4,1.25,1.25]};
simnames = {'11-disp-standard'};

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
        zUL = 4.0 %settings.zUL;
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
        
        %         gamma_13 =  1/Tdeph_1; %% dephsing of the resonant tunneling transition
        %         gamma_32 =  1/2*(G(ULL) + G(LLL)) + 1/Tdeph_2; % dephasing of the optical transision...
        %         gamma_12 =  1/2*(G(INJ) + G(LLL)) + 1/Tdeph_3;
        %
        
        dE13 = -1i*E13 - gamma_13; %
        dE32 = +1i*(E0 - E32) - gamma_32; %
        dE12 = +1i*(E0 - E12)- gamma_12; %
        
        % in the new system of units (mm/picosecond)
        
        gamma = 10^-11; % [m/V^2]
        
        alpha = gamma*(Constants('hbar')/( settings.zUL*10^-9*Constants('q0') ))^2;
        alpha = alpha*(10^24/(10^3));%renormalize :)
        
        n2 = 10^-13;%m^2/W
        epsilon = n2*E0/(c*1E-3);
        epsilon = epsilon*(Constants('hbar')/( settings.zUL*10^-9*Constants('q0') ))^2;
        epsilon = epsilon*(10^24/10^3);
        
        
        Ncarriers = settings.dN*(100^3)*settings.Ld/settings.Lp; % cm^-3 --> m^-3; carrier density
        Overlap = 0.8 ;%settings.Overlap;  % overlap factor -> dimensionless
        
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
        
        NRT = 0.8;
        tEnd = NRT*T_R; % end time in tps
        plotCtr = 1000; %% set on how many iterations should the program plot the intensity
        pulse_len = tEnd; % how many ps should we record the pulse for
        
        if(init > 0)
            
            x_0 = Ltot/7; tp =0.3;
            aE_in = @(z,time) exp(-(time-(z-x_0)/c).^2/tp^2);
            U = zeros(N,1);
            dw = 2*pi/T_R; sig = 2;
            %             nmods = 50;
            %             for m = -nmods/2:nmods/2-1
            %                 U = U+exp(1i*m*dw/c*(x-x_0));
            %             end
            U = aE_in(x,0);
            %normalize
            U = ampl*U/max(abs(U));
            
            U_t = zeros(5,1,1);  V_t = zeros(5,1,1);
            checkpoints = [1, round(N/5)+1, 2*round(N/5) , 3*round(N/5), 4*round(N/5)];
            ctr = 1;
            single_ctr = 1;
            idx = 1;
            
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
            
            
        end
        
        legendinfo = {'Intensity','inj','ull','lll'}; nlegend = length(legendinfo);
        for p = 1:length(idx_rest)
            legendinfo{nlegend+p} = [' lvl: ' num2str(idx_rest(p)) ];
        end
        
        probe_pt  = 5 ;
        
        while(t< tEnd)
            
            
            %%plot some of the results if neeed ariseth :D
            if(mod(iter_ctr,plotCtr) == 0)
                clc;
                display([settings.name ' simulation!']);
                display(['iteration Nr. = ' num2str(iter_ctr) ' @ RT = ' num2str(t/T_R)])
                plotyy(x,real(U),x,[r11,r33,r22,populations]);
                legend(legendinfo);
                title([ settings.name ' (MS + RNFD)  @ t = ' num2str(t)]);
                getframe;
            end
            
            rt = 1 ;
            U_t(1,rt,ctr) = U(checkpoints(1));
            U_t(2,rt,ctr) = U(checkpoints(2));
            U_t(3,rt,ctr) = U(checkpoints(3));
            U_t(4,rt,ctr) = U(checkpoints(4));
            U_t(5,rt,ctr) = U(checkpoints(5));
            ctr = ctr+1;
            
            
            %Update the round-trip counter!
            
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
            
            U = wave_solver.make_step(-1i*c*n32,-1i*c*n32_t,losses,dt);
            %set the boundaries... and obtain the final solution
            U = wave_solver.set_bdry(U(end),'no');
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
        
        
        mkdir(simnames{tau_idx});
        savename = [settings.name '_DA_' simnames{tau_idx} '_' strrep(num2str(ampl),'.','p') '_RING_' num2str(N) '_' strrep(num2str(Tdeph_1),'.','p') '_'  strrep(num2str(Tdeph_2),'.','p') '_'  strrep(num2str(Tdeph_3),'.','p')];
        save([simnames{tau_idx} '\' savename]);
        
    end
end
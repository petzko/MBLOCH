clear; clc;close all;
scenariofile = 'RTSLOWLIGHT_PROTOTYPE5-5e16doping.sim';
simfile = 'RTSLOWLIGHT_PROTOTYPE5-5e16doping.sim';

amplitudes = [0.001];
taus = {[1000,1000,1000]};% pure dephasing-lifetimes 
simnames = {'RTSLOWLIGHT5-SDA_5e16doping'};

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
        % hbar in eV-ps
        hbar = Constants('hbar',{'time',tch})/Constants('q0');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEVICE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%
        Ltot = 5; % mm
        n = settings.nTHz; 
        c=c_0; T_R = Ltot/c; f_R = 1/T_R;
        deg = 4.5;
        %slow light params
        transition_factor = +1;
        lossfactor =0;
        %optically thick medium
        dopingscale = +1;
        couplingfactor = 1 ;
        dephasing_factor = 1;
        %Courant Factor
        CF =1;
        %population inversion scheme
        J_rate = 0; %10 THz?  
        finite_support = false;
        
%%%%%%%%%%%%%%%% simuation parameters %%%%%%%%%%%%%%%%%
        
        %grid size in x direction
        N =2000; %nr of points in x direction
        x = linspace(0,Ltot,N)';
       
        dx = x(2) - x(1); dt = CF*dx/c;
       
        %number of regions 
        NR = 100; 
        % total length of slowling regions
        SL =  Ltot*0.3; 
        dW = SL/NR;
        
        %length of free space
        FL = (Ltot-SL);
        %width of a single Free space region 
        dF = FL/(NR+1); 
        
        %free-space slow-down free-space slow-down free-space
        core_i = zeros(length(x),NR);
        core = 0*x;
        for i = 1:NR
            core_i(:,i) = (x>= (i*dF+(i-1)*dW)  & x<= (i*dF+i*dW));
            core = core + core_i(:,i);
        end
        %interaction area ! 
        i_core = (x>(dF) & x<=(NR*dF+NR*dW))
        
         
        HTB = settings.HTB;
        NLVLS = sqrt(length(HTB)); % nr of levels to consider !
        HTB = reshape(HTB,NLVLS,NLVLS).';
        
        Spin_ = 3; Ex_ = 2; Grnd_ = 1;
        %there is no depopulation level and no ignore level
        IGNORELEVEL = -1; 
        DEPOP = -1;
        
        E_s = HTB(Spin_,Spin_)/hbar;
        E_e = HTB(Ex_,Ex_)/hbar;
        E_g = HTB(Grnd_,Grnd_)/hbar; % in rad/ps
        E_s = E_e ;
        %%%% energies: E1 > E3 > E2 Optical transition is 3->2 TUNNELING transition is 1->3
        
        
        O_se = couplingfactor*HTB(Spin_,Ex_)/hbar; % in rad/ps
        E_se = E_s-E_e; %rad/ps; 1->3 traisition freq
        E_sg = E_s-E_g; %rad/ps; 1->2 transition freq
        E_eg = E_e-E_g; %rad/ps  3->2 transition freq (optical transition)
    
        % central frequency and wave number!
        E0 = E_eg; % central OPTICAL frequency.
        
        %gain recovery times
        W_absorb = settings.Wmtx;
        NLVLS = sqrt(length(W_absorb)); % nr of levels to consider !
        %reshape into a matrix
        W_absorb = reshape(W_absorb,NLVLS,NLVLS).';
        for j = 1:NLVLS
            W_absorb(j,j) = 0;
        end
        
        W_amplify = W_absorb;
        W_amplify(1,2) = 1e-2;
        W_amplify(2,1) = 1e-5;
        W_amplify(2,3) = 1e-5;
        W_amplify(3,2) = 1e-5;
        
        temp = W_amplify(2,:); 
%         W_amplify(2,:) = W_amplify(1,:); 
%         W_amplify(1,:) = temp;
        W = zeros(NLVLS,NLVLS,2); 
        W(:,:,1) = W_absorb; 
        W(:,:,2) = W_amplify;
        
        %extract the indices of the rest of the laser levels, which will be
        %included in the system as simple rate equations
        idx_rest = setdiff(1:NLVLS,[Spin_,Ex_,Grnd_,DEPOP,IGNORELEVEL]);
        N_rest = length(idx_rest); %take the number of residual levels!
        
        G_ = zeros(NLVLS,2);
        zeroFCTR = 1;
        
        %absorption layer lifetimes
        G_(Spin_,1) = sum(W(Spin_,setdiff(1:NLVLS,[Spin_,DEPOP,IGNORELEVEL]),1));
        for lvl = setdiff(1:NLVLS,[Spin_,DEPOP,IGNORELEVEL])
            G_(lvl,1) = sum(W(lvl,setdiff(1:NLVLS,[IGNORELEVEL]),1));
        end
        
        %amplification layer lifetimes
        G_(Spin_,2) = sum(W(Spin_,setdiff(1:NLVLS,[Spin_,DEPOP,IGNORELEVEL]),2));
        for lvl = setdiff(1:NLVLS,[Spin_,DEPOP,IGNORELEVEL])
            G_(lvl,2) = sum(W(lvl,setdiff(1:NLVLS,[IGNORELEVEL]),2));
        end
      
        G = zeros(length(core),NLVLS); 
        
        for j = 1:NLVLS
            G(:,j) = G_(j,1).*core + G_(j,2).*(1-core);
        end
        
        Tdeph_se = taus{tau_idx}(1);
        Tdeph_eg = taus{tau_idx}(2);
        Tdeph_sg = taus{tau_idx}(3);
        
        % % % % Added Pure Dephasing
        gamma_se = dephasing_factor*(  1/2*(G(:,Spin_) + G(:,Ex_)) + 1/Tdeph_se); %% dephsing of the resonant tunneling transition
        gamma_eg = dephasing_factor*(  1/2*(G(:,Ex_) + G(:,Grnd_)) + 1/Tdeph_eg); % dephasing of the optical transision...
        gamma_sg = dephasing_factor*(  1/2*(G(:,Spin_) + G(:,Grnd_))+ 1/Tdeph_sg); % dephasing of the latest transition
        
%         gamma_sg = 0.*core + gamma_sg.*(1-core);
        
        dE_se = -1i*E_se - gamma_se; %
        dE_eg = +1i*(E0 - E_eg) - gamma_eg; %
        dE_sg = +1i*(E0 - E_sg)- gamma_sg; %
        
        Ncarriers = dopingscale*settings.dN*(100^3)*settings.Ld/settings.Lp; % cm^-3 --> m^-3; carrier density
        Overlap = 0.95;  % overlap factor -> dimensionless
        %the plasma freq squared ;)
        Wp= 1e-12*Ncarriers*Overlap*(deg*1E-9*Constants('q0'))^2/(Constants('eps0')*Constants('hbar')*n^2);
        
        display(['Plasma freq (THz):' num2str(E0*Wp) ''] );        
        display(['Coupling strength ^2 (THz):' num2str(O_se^2)])
        
        %calculate the normalization constant, i.e. the number on which we divide
        %all quantities ( except the electric field envelope) to simplify
        %our initial system. it is important as it determines the initial value of
        %the overall electron population inside the system-> a quantity that shall
        %be perserved throughout the whole simulaiton ! !
        trace_rho = ((E0*1E12*Ncarriers*Overlap*((deg*1E-9*Constants('q0'))^2))/(Constants('eps0')*n*Constants('c')*Constants('hbar')))/(1/(lch*tch));
        
        %cavity loss l_0 in (cm^-1) --> l_0*100 in (m^-1) --> 1 mm^-1
        l_0 = settings.loss*100/(1/lch);
        
        %
        
        k0 = E0/c;
        
        %%%% specify some of the mainloop control parameters %%%%
        iter_per_rt = round(T_R/dt);
        
        NRT = 500;
        tEnd = NRT*T_R; % end time in tps
        plotCtr = 1000; %% set on how many iterations should the program plot the intensity
        pulse_len = tEnd; % how many ps should we record the pulse for
        
        if(init > 0)
            
            FWHM = 30*2/sqrt(2); %psec
            t0 = +50;
            sig = FWHM/(2*sqrt(2*log(2))); % std. dev
            
            gauss = @(t,mu,s) ampl*exp(-(t-mu).^2/(2*s^2));            
            
            NBits =  8; 
            NUM = 125; 
            code = de2bi(NUM,NBits); 
            message = @(t) 0;
            
            for i = 1:NBits
                message = @(t) message(t) + code(i)*gauss(t,t0+i*3*FWHM,sig)
            end
            
            bcar = @(t,t_i,t_f) ampl*(t>t_i & t< t_f);
            
            if finite_support
                pulsefunc = @(t,mu) bcar(t,mu,FWHM+mu);
            else
                pulsefunc = @(t,mu) gauss(t,mu,sig); 
            end
            
            U = zeros(N,1);
            
            U_in =  zeros(round(tEnd/dt),1);
            U_out =  zeros(round(tEnd/dt),1);
            
            ctr = 1;
            single_ctr = 1;
            idx = 1;
            
            trace_rho_amp = 5*trace_rho;
            % population vectors together with first derivative vectors
            r_ee = trace_rho*0*ones(N,1); r_ee(1) = r_ee(end);
            r_gg = trace_rho*1*ones(N,1).*core +  trace_rho_amp*ones(N,1).*(1-core); r_gg(1) = r_gg(end);
            r_ss = trace_rho*0*ones(N,1); r_ss(1)=r_ss(end);
            
            %rate populations
            populations = zeros(N,N_rest);
            
            % coherence terms
            r_se = 1E-20*((rand(N,1)-0.5) +1i*(rand(N,1)-0.5)); r_se(1) = r_se(end);
            n_eg = 1E-20*((rand(N,1)-0.5) +1i*(rand(N,1)-0.5)); n_eg(1) = n_eg(end);
            n_sg = 1E-20*((rand(N,1)-0.5) +1i*(rand(N,1)-0.5)); n_sg(1) = n_sg(end);
            
            t = 0; nr_steps = 5;
            
            iter_ctr = 0;
            iterperrecord = 1 ; % nr of iterations after which to record new statistics!
            
            r_ss_solver = MS(nr_steps,N,[],r_ss); r_gg_solver = MS(nr_steps,N,[],r_gg);
            r_ee_solver = MS(nr_steps,N,[],r_ee); n_sg_solver = MS(nr_steps,N,[],n_sg);
            r_se_solver = MS(nr_steps,N,[],r_se); n_eg_solver = MS(nr_steps,N,[],n_eg);
            
            rate_eqn_solvers = [];
            for p=1:N_rest
                rate_eqn_solvers{p} =  MS(nr_steps,N,[],populations(:,p));
            end
            wave_solver = RNFDSolver(N,dx,c/abs(c),abs(c), U);
%             wave_solver = LaxSolver(N,dx,c/abs(c),abs(c), U);
        end
      
        
        probe_pt  = 5 ;
        T_decay = 50000000; 
        tms = 1;
        frame_ctr = 1; 
        while(t< tEnd)

            
            %%plot some of the results if neeed ariseth :D
            if(mod(iter_ctr,100) == 0)
                clc;
                r = mean(r_ee - r_gg);
                display(['(n-1)/n*(Coupling strength ^2 /w0)(THz):' num2str((n-1)/n*O_se^2/E0)])
                display(['r x Plasma Freq (THz):' num2str(r*Wp) ''] );
                display(['(Coupling strength ^2 /w0)(THz):' num2str(O_se^2/E0)])
                display([settings.name ' simulation!']);
                display(['iteration Nr. = ' num2str(iter_ctr) ' @ RT = ' num2str(t/T_R)])
                Op = Ncarriers*Overlap*(Constants('q0')*deg*1e-9)^2/(n^2*Constants('eps0')*Constants('hbar'));
%                 r_inv = mean(r_ee-r_gg)/trace_rho;
%                 ng1 = n*(1-r_inv*E0*1e12*Op*((O_se*1e12)^2-(gamma_sg*1e12)^2)/((O_se*1e12)^2+(gamma_sg*1e12)*(gamma_eg*1e12))^2);
%                 ng2 = n*(1-0.5*r_inv*E0*1e12*Op*((O_se*1e12)^2-(gamma_sg*1e12)^2)/((O_se*1e12)^2+(gamma_sg*1e12)*(gamma_eg*1e12))^2);
%                 
%                 display(['ng1 = ' num2str(ng1) '; ng2 = ' num2str(ng2) '; nthz = ' num2str(n)])
%                 display(['Vg1 = ' num2str(c_0/ng1) '; Vg2 = ' num2str(c_0/ng2) '; Vph = ' num2str(c_0/n)])
                subplot(2,1,1);
                ax = plot(x,abs(U).^2,'Color','r');
%                 xlim([dF,dF+dW])
                hold on;
                fobj = fill(x,core*ampl^2,'b'); 
                set(fobj,'EdgeColor','none','FaceAlpha',0.4);
                hold off;
                subplot(2,1,2);
                ax = plot(x,[core,(r_ee-r_gg).*core/trace_rho,(r_ee-r_gg).*(1-core)/trace_rho_amp]);
%               title([ settings.name ' (MS + RNFD)  @ t = ' num2str(t)]);
                getframe;
%                 if (mod(iter_ctr,100) == 0)
%                     frames(frame_ctr) = getframe; 
%                     frame_ctr = frame_ctr +1; 
% %                     pause;
%                     
%                 else
%                     getframe;
%                 end
                 
            end
            
            rt = 1;
            U_in(ctr) = U(1);
            U_out(ctr) = U(end);
            ctr = ctr+1;
            
            
            %Update the round-trip counter!
            
            %%%%%% BLOCH PART %%%%%%
            %%%% POPULATIONS
            rss_t = 1i*O_se*(r_se-conj(r_se)) +...
                     W(Ex_,Spin_,1)*r_ee.*core +...
                     W(Ex_,Spin_,2)*r_ee.*(1-core) +...
                    (W(Grnd_,Spin_,1) +J_rate)*r_gg.*core +...
                    (W(Grnd_,Spin_,2) +J_rate)*r_gg.*(1-core) -...
                    G(:,Spin_).*r_ss;
            for p = 1:N_rest
                p_glob_idx = idx_rest(p);
                rss_t = rss_t + (W(p_glob_idx,Spin_,1).*core + W(p_glob_idx,Spin_,2).*(1-core)).*populations(:,p);
            end
            r_ss_solver.make_step(rss_t,dt);
            
            lmInteraction = conj(U).*n_eg;
            ree_t = -1i*O_se*(r_se-conj(r_se)) + ...
                    1i/2.*(lmInteraction-conj(lmInteraction)) + ...
                    core.*r_ss*(W(Spin_,Ex_,1)) + ...
                    (1-core).*r_ss*(W(Spin_,Ex_,2)) + ...
                    core.*r_gg*W(Grnd_,Ex_,1) +...
                    (1-core).*r_gg*W(Grnd_,Ex_,2) - G(:,Ex_).*r_ee;
            for p = 1:N_rest
                p_glob_idx = idx_rest(p);
                ree_t = ree_t + (W(p_glob_idx,Ex_,1).*core+W(p_glob_idx,Ex_,2).*(1-core)).*populations(:,p);
            end
            r_ee_solver.make_step(ree_t,dt);
            
            rgg_t = -1i/2.*(lmInteraction-conj(lmInteraction)) + ...
                core.*r_ss*W(Spin_,Grnd_,1) + ...
                (1-core).*r_ss*W(Spin_,Grnd_,2)+ ...
                core.*r_ee*W(Ex_,Grnd_,1) + ...
                (1-core).*r_ee*W(Ex_,Grnd_,2) - ...
                (G(:,Grnd_)+J_rate).*r_gg ;
           
            for p = 1:N_rest
                p_glob_idx = idx_rest(p);
                rgg_t = rgg_t + (W(p_glob_idx,Grnd_,1).core+ ...
                    W(p_glob_idx,Grnd_,2).(1-core))*populations(:,p);
            end
            r_gg_solver.make_step(rgg_t,dt);
            
            for p = 1:N_rest
                p_glob_idx =idx_rest(p);
                rpp_0_t = W(Spin_,p_glob_idx,1)*r_ss.*core + ... 
                          W(Spin_,p_glob_idx,2)*r_ss.*(1-core) + ...
                          W(Ex_,p_glob_idx,1)*r_ee.*core + ... 
                          W(Ex_,p_glob_idx,2)*r_ee.*(1-core) + ... 
                          W(Grnd_,p_glob_idx,1)*r_gg.*core + ...
                          W(Grnd_,p_glob_idx,2)*r_gg.*(1-core) + ...
                          - G(:,p_glob_idx).*populations(:,p);
                
                %add the rate equations part!
                for j =1:N_rest % nr
                    if j~= p
                        j_glob_idx =idx_rest(j);
                        rpp_0_t = rpp_0_t + (W(j_glob_idx,p_glob_idx,1).*core + ...
                        W(j_glob_idx,p_glob_idx,2).*(1-core)).*populations(:,j);
                    end
                end
                rate_eqn_solvers{p}.make_step(rpp_0_t,dt);
                
            end
            
            %%% coherences! r13 n32 n12
            rse_t = dE_se.*r_se + 1i*O_se*(r_ss - r_ee) + 1i/2*conj(U).*n_sg;
            r_se_solver.make_step(rse_t,dt);
            
            neg_t = dE_eg.*n_eg + 1i/2*U.*(r_ee-r_gg) - 1i*O_se*n_sg;
            n_eg_solver.make_step(neg_t,dt);
            
            n12_t = dE_sg.*n_sg + 1i/2*U.*r_se - 1i*O_se*n_eg;
            n_sg_solver.make_step(n12_t,dt);
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            losses = -c*l_0*lossfactor*ones(N,1);
         
            U = wave_solver.make_step(-1i*c*n_eg.*i_core*transition_factor,-1i*c*neg_t.*i_core*transition_factor,losses,dt);
%              U = wave_solver.make_step(-1i*c*n_eg*0.*i_core*transition_factor,-1i*c*n_eg*0.*i_core*transition_factor,losses,dt);
            %set the boundaries... and obtain the final solution
            %set the boundaries... and obtain the final solution
            U = wave_solver.set_bdry(message(t),'no');
            %%%%%%%%%%%%%%%%%%%%%%%%
            r_ss = r_ss_solver.get_latest_solution();
            r_ee = r_ee_solver.get_latest_solution();
            r_gg = r_gg_solver.get_latest_solution();

%             r_ss = 0*core + (1-core).*r_ss_solver.get_latest_solution();
%             r_ee = 0*core + (1-core).*r_ee_solver.get_latest_solution();
%             r_gg = trace_rho.*core+(1-core).*r_gg_solver.get_latest_solution();
            
            for p = 1:length(idx_rest)
                populations(:,p) = rate_eqn_solvers{p}.get_latest_solution();
            end
            
            r_se = r_se_solver.get_latest_solution();
            n_eg = n_eg_solver.get_latest_solution();
            n_sg = n_sg_solver.get_latest_solution();
            
            t = t+dt;
            iter_ctr = iter_ctr + 1;
        end
        
        
        mkdir(simnames{1});
        savename = [simnames{1} 'DA_version' num2str(tau_idx)];
        save([simnames{1} '\' savename]);
        
    end
end
function [ savename ] = runsim_MB_TL_mlvl( scenarioFile,simDataFile,ratesfile,varargin )


close all;

lvar = length(varargin);

%handle user input!
if(length(varargin)>0)
    if mod(lvar,2) ~=0
        display('incorrect number of input arguments. Aborting!');
        return;
    end
    for i = 1:2:lvar
        name = varargin{i}; val  = varargin{i+1};
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

%parse all input files and load the scatterin rates file !
settings = parseSimParams(scenarioFile);
settings = parseSimData(simDataFile,settings);
load(ratesfile);

tch = settings.tch; % characteristic time..
lch = settings.lch; % characteristic length...

%speed of light in vacuum: 299792458m/s --> in tch/lch
c_0 = Constants('c',{'time',tch},{'length',lch});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEVICE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%
%cavity length  -->
Ltot = settings.Ltot; % mm

%phase velocity inside the medium ( in mm per picosecond ... )
n = settings.n;  c = c_0/n; T_R = 2*Ltot/c; f_R = 1/T_R;

%%%%dipole mtx elements (in Cnm)
zUL = settings.zUL;
% hbar in eV-ps
hbar = Constants('hbar',{'time',tch})/Constants('q0');

if abs(HTB(1,3))>abs(HTB(2,3))
    INJ = 1; IGNORELEVEL = 7; DEPOP = 6;
else
    INJ =2; IGNORELEVEL = 6; DEPOP = 7;
end
ULL = 3; LLL = 4;

if (settings.shb > 0)
    shb = 1;
else
    shb = -1;
end

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
    x_0 = Ltot/7; tp =1;
    aE_in = @(z,time) exp(-(time-(z-x_0)/c).^2/tp^2);
    
    U = aE_in(x,0);
    %normalize
    U = 3*U/max(abs(U));
    V = 0*U;
    
    
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
    V_TL_t =1;
    J_TL_t =1;
    
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
    
    
    %%% transmission line params
    % transmission line params:
    bias = 11/1E1; %initial bias in kV/mm;
    mag = bias*1E-1;
    v_0 =@(tm) bias + mag*sin(2*pi*1/T_R*tm);
    v_TL_new = v_0(0)*ones(N,1); % transmission line voltage per unit length (in units kV/mm);
    
    
    NLVLS = 7; %nr of levels to consider !
    
    idx_rest = setdiff(1:NLVLS,[INJ,ULL,LLL,DEPOP,IGNORELEVEL]);
    N_rest = length(idx_rest); %take the number of residual levels!
    G = zeros(settings.N,NLVLS);
    v_TL_old = v_TL_new;
    
    J_TL = zeros(N,1); % QCL curr density (in units A/mm^2)
    rates = zeros(N,1); rates_t = zeros(N,1);
    %%% curr density calculation:
    rates = (r110).*W_fit{INJ,DEPOP}(v_TL_new) + (r220).*W_fit{LLL,DEPOP}(v_TL_new) + (r330).*W_fit{ULL,DEPOP}(v_TL_new);
    for p = 1:N_rest
        p_glob_idx = idx_rest(p);
        rates = rates + (populations(:,p)).*W_fit{p_glob_idx,DEPOP}(v_TL_new);
    end
    
    %convert the curr density in A/m^2 into A/mm^2 by dividing by 1E6;
    J_TL = (Constants('q0')*(settings.Lp*1E-9)*Ncarriers/trace_rho*1E12*rates)/1E6; %in A/mm^2
    v_tmp = v_TL_new;
    eps_ch = 1E15*Constants('eps0');    mu_ch  = 1E3*Constants('mu0');
    %%%
    
end

iterperrecord = 1;
recordingiter  = round(recordingduration/iterperrecord/dt);
padsize = recordingiter-length(E_p);

%preallocate memory for optical field, curr density voltage wave, etc.
%strorage...
E_p= padarray(E_p,padsize,'post'); E_m = padarray(E_m,padsize,'post');

V_TL_t = padarray(V_TL_t,padsize,'post'); J_TL_t = padarray(J_TL_t,padsize,'post');

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


info.settings = settings;
info.cavity = 'FP';
info.Ltot = Ltot;
info.N = N;
info.SIMTYPE = 'WITHOUT DISP. COMP';
info.l_0 = l_0;



while(t< tEnd)
    
    %%% begin setting up TL params!
    for el = setdiff(1:NLVLS,[INJ,DEPOP,IGNORELEVEL])
        G(:,INJ) = G(:,INJ)+W_fit{INJ,el}(v_TL_new);
    end
    
    for lvl = setdiff(1:NLVLS,[INJ,DEPOP,IGNORELEVEL])
        for el = 1:NLVLS
            G(:,lvl) = G(:,lvl) +  W_fit{lvl,el}(v_TL_new);
        end
    end
    
    %%%%% correctly setup the energies, the anticrossings and the scattering rates...
    %obtain the energies of the core levels for the simulation
    E1 = E_fit{INJ}(v_TL_new)/hbar;
    E3 = E_fit{ULL}(v_TL_new)/hbar;
    E2 = E_fit{LLL}(v_TL_new)/hbar; % in rad/ps
    %%%% energies: E1 > E3 > E2 Optical transition is 3->2 TUNNELING transition
    %%%% is 1->3
    O13 = O_13_fit(v_TL_new)/hbar; % in rad/ps
    E13 = E1-E3; %rad/ps; 1->3 traisition freq
    E12 = E1-E2; %rad/ps; 1->2 transition freq
    E32 = E3-E2; %rad/ps  3->2 transition freq (optical transition)
    
    %extract the indices of the rest of the laser levels, which will be
    %included in the system as simple rate equations
    
    % % % % Added Pure Dephasing
    if(settings.deph>0)
        gamma_13 = 0.5*(G(:,INJ)+G(:,ULL))+1/settings.Tdeph_1; %% dephsing of the resonant tunneling transition
        gamma_32 = 0.5*(G(:,LLL)+G(:,ULL))+1/settings.Tdeph_2; % dephasing of the optical transision...
        gamma_12 = 0.5*(G(:,LLL)+G(:,INJ))+1/settings.Tdeph_3; % dephasing of the latest transition
    else
        gamma_13 = 0.5*(G(:,INJ)+G(:,ULL)); %% dephsing of the resonant tunneling transition
        gamma_32 = 0.5*(G(:,LLL)+G(:,ULL)); % dephasing of the optical transision...
        gamma_12 = 0.5*(G(:,LLL)+G(:,INJ)); % dephasing of the latest transition
    end
    
    dE13 = -1i*E13 - gamma_13; %
    dE32 = +1i*(E0 - E32) - gamma_32; %
    dE12 = +1i*(E0 - E12)- gamma_12; %
        
    %%%%% end of setting up the TL PARAMS %%%%%
    
    
    if(mod(iter_ctr+1,checkptIter) == 0 )
        checkptname = ['CHCKPT_' settings.name '_' settings.scenario '_N_' num2str(settings.N) '_FP'];
        save(checkptname);
    end
    
    %%plot some of the results if neeed ariseth :D
    if(mod(iter_ctr,plotCtr) == 0)
        clc;
        info.iter_ctr = iter_ctr;
        info.RT = t/T_R;
        intensity = U.*conj(U) + V.*conj(V) ;
        info.maxInt  =  max(intensity);
        printINFO(info);
        
        plot(x,[abs(U),abs(V)]);
        getframe;
    end
    %%%% obtain the field, field intensity and the total population at position "idx" ...
    
    if((t >= tEnd - recordingduration) && mod(iter_ctr,iterperrecord) == 0)
        
        %store fields info
        E_p(ctr)= U(idx);  E_m(ctr)= V(idx);
        V_TL_t(ctr) = v_TL_new(idx);
        J_TL_t(ctr) = J_TL(idx);
        
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
    r110_t = 1i*O13.*(r130-conj(r130)) +((W_fit{ULL,INJ}(v_TL_new) + W_fit{ULL,DEPOP}(v_TL_new)).*r330 + ( W_fit{LLL,INJ}(v_TL_new)+W_fit{LLL,DEPOP}(v_TL_new) ).*r220 - G(:,INJ).*r110;
    for p = 1:N_rest
        p_glob_idx = idx_rest(p);
        r110_t = r110_t + (W_fit{p_glob_idx,INJ}(v_TL_new)+W_fit{p_glob_idx,DEPOP}(v_TL_new)).*populations(:,p);
    end
    rates_t = r110_t.*W_fit{INJ,DEPOP}(v_TL_new);
    r110_solver.make_step(r110_t,dt);
    
    lmInteraction = conj(U).*n32p + conj(V).*n32m;
    r330_t = 1i*O13.*(conj(r130) - r130) +1i/2.*(lmInteraction-conj(lmInteraction)) +  r110.*W_fit{INJ,ULL}(v_TL_new) + r220.*W_fit{LLL,ULL}(v_TL_new) - G(:,ULL).*r330;
    for p = 1:N_rest
        p_glob_idx = idx_rest(p);
        r330_t = r330_t + W_fit{p_glob_idx,ULL}(v_TL_new).*populations(:,p);
    end
    rates_t = rates_t + r330_t.*W_fit{ULL,DEPOP}(v_TL_new);
    r330_solver.make_step(r330_t,dt);
    
    r220_t = -1i/2.*(lmInteraction-conj(lmInteraction)) + r110.*W_fit{INJ,LLL}(v_TL_new) + r330.*W_fit{ULL,LLL}(v_TL_new) - G(:,LLL).*r220;
    for p = 1:N_rest
        p_glob_idx = idx_rest(p);
        r220_t = r220_t + W_fit{p_glob_idx,LLL}(v_TL_new).*populations(:,p);
    end
    rates_t = rates_t + r220_t.*W_fit{LLL,DEPOP}(v_TL_new);
    r220_solver.make_step(r220_t,dt);
    
     for p = 1:N_rest
        p_glob_idx =idx_rest(p);
        rpp_0_t = W_fit{INJ,p_glob_idx}(v_TL_new).*r110+W_fit{ULL,p_glob_idx}(v_TL_new).*r330+W_fit{LLL,p_glob_idx}(v_TL_new).*r220 - G(:,p_glob_idx).*populations(:,p);
        %add the rate equations part!
        for j =1:N_rest % nr
            if j~= p
                j_glob_idx =idx_rest(j);
                rpp_0_t = rpp_0_t + W_fit{j_glob_idx,p_glob_idx}(v_TL_new).*populations(:,j);
            end
        end
        rates_t = rates_t + rpp_0_t.*W_fit{p_glob_idx,DEPOP}(v_TL_new);
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
        r11p_t = 1i*O13.*(r13p-conj(r13m)) + (W_fit{ULL,INJ}(v_TL_new)+W_fit{ULL,DEPOP}(v_TL_new)).*r33p+ (W_fit{LLL,INJ}(v_TL_new)+W_fit{LLL,DEPOP}(v_TL_new)).*r22p - (G(:,INJ)+diffusion).*r11p;
        r11p_solver.make_step(r11p_t,dt);
        
        %%% r33+
        r33p_t = 1i*O13.*(conj(r13m)-r13p)+1i/2*(conj(V).*(n32p) -(U).*conj(n32m)) +  W_fit{INJ,ULL}(v_TL_new).*r11p + W_fit{LLL,ULL}(v_TL_new).*r22p - (G(:,ULL)+diffusion).*r33p;
        r33p_solver.make_step(r33p_t,dt);
        
        %%% r22+
        r22p_t = -1i/2*(conj(V).*(n32p) -(U).*conj(n32m)) + W_fit{INJ,LLL}(v_TL_new).*r11p + W_fit{ULL,LLL}(v_TL_new).*r33p - (G(:,LLL)+diffusion).*r22p;
        r22p_solver.make_step(r22p_t,dt);
        
        %%% r13+
        r13p_t = (dE13-diffusion)*r13p + 1i*O13*(r11p-r33p) +1i/2*(conj(V).*n12p);
        r13p_solver.make_step(r13p_t,dt);
        
        %%% r13-
        r13m_t = (dE13-diffusion)*r13m + 1i*O13*conj(r11p-r33p) + 1i/2*(conj(U).*n12m);
        r13m_solver.make_step(r13m_t,dt);
    end
    
     %%%%% Begin TRANSMISSION LINE EQUATIONS %%%%%%%%%%%%%
     
    %%% curr density derivative calculation:
    rates = (r110).*W_fit{INJ,DEPOP}(v_TL_new) + (r220).*W_fit{LLL,DEPOP}(v_TL_new) + (r330).*W_fit{ULL,DEPOP}(v_TL_new);
    for p = 1:N_rest
        p_glob_idx = idx_rest(p);
        rates = rates + (populations(:,p)).*W_fit{p_glob_idx,DEPOP}(v_TL_new);
    end
    %convert the curr density in A/m^2 into A/mm^2 by dividing by 1E6;
    J_TL = (Constants('q0')*(settings.Lp*1E-9)*Ncarriers/trace_rho*1E12*rates)/1E6; %in A/mm^2
    J_TL_t = (Constants('q0')*(settings.Lp*1E-9)*Ncarriers/trace_rho*1E12*rates_t)/1E6;
    
    v_tmp(2:end-1) = v_TL_new(3:end)+v_TL_new(1:end-2)-v_TL_old(2:end-1)-dt^2/eps_ch/n^2*J_TL_t(2:end-1);
    
    %set bdry conditions for v at right end! 
    v_tmp(end) = v_tmp(end-1); v_tmp(1) = v_0(t); 
    v_TL_old = v_TL_new; v_TL_new = v_tmp;  
    
    %%%%% END TRANSMISSION LINE EQUATIONS %%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    losses = -c*l_0.*ones(N,1);
    
    U = U_solver.make_step(-1i*c*n32p,-1i*c*n32p_t,losses,dt);
    V = V_solver.make_step(-1i*c*n32m,-1i*c*n32m_t,losses,dt);
    
    %set the boundaries... and obtain the final solution
    U = U_solver.set_bdry(V(1),'no');
    V = V_solver.set_bdry('no',U(N));
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    r110 = r110_solver.get_latest_solution();
    r330 = r330_solver.get_latest_solution();
    r220 = r220_solver.get_latest_solution();
    
    for p = 1:length(idx_rest)
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
savename = [settings.name '_' settings.scenario '_N_' num2str(settings.N) '_FP_' num2str(settings.simRT) ];
save(savename);
end


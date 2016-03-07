function [ savename ] = runsim_RING_RNFD_multilevelHTB( scenarioFile,simDataFile,varargin )

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

%parse all input files
settings = parseSimParams(scenarioFile);
settings = parseSimData(simDataFile,settings);

tch = settings.tch; % characteristic time..
lch = settings.lch; % characteristic length...

%speed of light in vacuum: 299792458m/s --> in tch/lch
c_0 = Constants('c',{'time',tch},{'length',lch});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEVICE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%
%cavity length  -->
Ltot = settings.Ltot; % mm
%phase velocity inside the medium ( in mm per picosecond ... )
n = settings.n;  c = c_0/n; T_R = Ltot/c; f_R = 1/T_R;

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
global_idx_rest = setdiff(1:NLVLS,[INJ,ULL,LLL,DEPOP,IGNORELEVEL]);
N_rest = length(global_idx_rest); %take the number of residual levels!

G = zeros(NLVLS,1);
% W(DEPOP,:) = zeros(NLVLS,1); W(INJ,DEPOP) = 0;
G(INJ) = sum(W(INJ,setdiff(1:NLVLS,[INJ,DEPOP,IGNORELEVEL])));
for lvl = setdiff(1:NLVLS,[INJ,DEPOP,IGNORELEVEL])
      G(lvl) = sum(W(lvl,setdiff(1:NLVLS,[IGNORELEVEL])));
end


Tdeph_1 = settings.Tdeph_1;
Tdeph_2 = settings.Tdeph_2;
Tdeph_3 = settings.Tdeph_3;

% % % % Added Pure Dephasing
if(settings.deph>0)
    gamma_13 = 0.5*(G(INJ)+G(ULL))+1/Tdeph_1; %% dephsing of the resonant tunneling transition
    gamma_32 = 0.5*(G(LLL)+G(ULL))+1/Tdeph_2; % dephasing of the optical transision...
    gamma_12 = 0.5*(G(LLL)+G(INJ))+1/Tdeph_3; % dephasing of the latest transition
else
    gamma_13 = 0.5*(G(INJ)+G(ULL)); %% dephsing of the resonant tunneling transition
    gamma_32 = 0.5*(G(LLL)+G(ULL)); % dephasing of the optical transision...
    gamma_12 = 0.5*(G(LLL)+G(INJ)); % dephasing of the latest transition
end

dE13 = -1i*E13 - gamma_13; %
dE32 = +1i*(E0 - E32) - gamma_32; %
dE12 = +1i*(E0 - E12)- gamma_12; %


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

tEnd = settings.simRT*T_R; % end time in tps
plotCtr = settings.plotCtr; %% set on how many iterations should the program plot the intensity
recordingduration = settings.recordRT*T_R; % how many ps should we record the pulse for

if(init > 0)

       
    clc; t = 0;
    U = zeros(N,1); 
    
    % population vectors together with first derivative vectors
    r11 = trace_rho*(ones(N,1)*1/3); r11(1) = r11(end);
    r33 = trace_rho*(ones(N,1)*1/3); r33(1) = r33(end);
    r22 = trace_rho*(ones(N,1)*1/3); r22(1)=r22(end);
    
    %rate populations
    populations = zeros(N,N_rest);
    
    % coherence terms together with first derivative vectors!
    r13 =0*((ones(N,1)-0.5) +1i*(ones(N,1)-0.5)); r13(1) = r13(end);
    n32 =  1E-15*((rand(N,1)-0.5) +1i*(rand(N,1)-0.5)); n32(1) = n32(end);
    n12 =  1E-15*((rand(N,1)-0.5) +1i*(rand(N,1)-0.5)); n12(1) = n12(end);
    
    t = dt; nr_steps = settings.nr_steps;
    
    idx = 1; % the index of the predefined point we sample the feild at! 
    iter_ctr = 0;    ctr = 1;   iterperrecord = 1 ; % nr of iterations after which to record new statistics!
    
     %simulation info storage arrays -> preallocate 
    recordingiter  = round(recordingduration/iterperrecord/dt);
    E_p= zeros(1,recordingiter);  
    %store population info 
    r11_time = E_p;    r33_time = E_p;
    r22_time = E_p;     
    pop_time = zeros(N_rest,recordingiter);
    r13_time = E_p; n32_time = E_p; n12_time =E_p;

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
for p = 1:length(global_idx_rest)
    legendinfo{nlegend+p} = [' lvl: ' num2str(global_idx_rest(p)) ];
end

info.settings = settings; 
info.cavity = 'RING';
info.Ltot = Ltot; 
info.N = N;
info.l_0 = l_0; 
info.SIMTYPE = 'WITHOUT DISPERSION COMPENSATION';


while(t< tEnd)
    
    
    %%plot some of the results if neeed ariseth :D
    if(mod(iter_ctr,plotCtr) == 0)
        clc;
        info.iter_ctr = iter_ctr; 
        info.RT = t/T_R;
        intensity = U.*conj(U) ;
        info.maxInt  =  max(intensity);
        printINFO(info); 
%         plotyy(x,real(U),x,[r11,r33,r22,populations]);
%         legend(legendinfo);
%         title([ settings.name ' (MS + RNFD)  @ t = ' num2str(t)]);
%         getframe;
    end
    
     if(t >= tEnd - recordingduration && mod(iter_ctr,iterperrecord) == 0)
        E_p(ctr)= U(idx);   
        %store population info 
        r11_time(ctr)= r11(idx);
        r33_time(ctr)= r33(idx);
        r22_time(ctr) = r22(idx);
        %calculate total populations
        for p = 1:N_rest
            pop_time(p,ctr) = populations(idx,p);
        end
        
        r13_time(ctr) = r13(idx); n32_time(ctr) = n32(idx); 
        n12_time(ctr) = n12(idx); 
        
        ctr = ctr+1;
    end
    
    %%%%%% BLOCH PART %%%%%%
    %%%% POPULATIONS
    r11_t = 1i*O13*(r13-conj(r13)) +( W(ULL,INJ) + W(ULL,DEPOP) )*r33 + ( W(LLL,INJ)+W(LLL,DEPOP) )*r22 - G(INJ)*r11;
    for p = 1:N_rest
        p_glob_idx = global_idx_rest(p); 
        r11_t = r11_t + (W(p_glob_idx,INJ)+W(p_glob_idx,DEPOP))*populations(:,p);
    end
    r11_solver.make_step(r11_t,dt);
    
    lmInteraction = conj(U).*n32;
    r33_t = -1i*O13*(r13-conj(r13)) +1i/2.*(lmInteraction-conj(lmInteraction)) + r11*W(INJ,ULL) + r22*W(LLL,ULL) - G(ULL)*r33;
    for p = 1:N_rest
        p_glob_idx = global_idx_rest(p); 
        r33_t = r33_t + W(p_glob_idx,ULL)*populations(:,p);
    end
    r33_solver.make_step(r33_t,dt);
    
    r22_t = -1i/2.*(lmInteraction-conj(lmInteraction)) + r11*W(INJ,LLL) + r33*W(ULL,LLL) - G(LLL)*r22;
    for p = 1:N_rest
        p_glob_idx = global_idx_rest(p); 
        r22_t = r22_t + W(p_glob_idx,LLL)*populations(:,p);
    end
    r22_solver.make_step(r22_t,dt);
    
    for p = 1:N_rest
        p_glob_idx =global_idx_rest(p); 
        rpp_0_t = W(INJ,p_glob_idx)*r11+W(ULL,p_glob_idx)*r33+W(LLL,p_glob_idx)*r22 - G(p_glob_idx)*populations(:,p);
        
        %add the rate equations part!        
        for j =1:N_rest % nr 
            if j~= p
                j_glob_idx =global_idx_rest(j); 
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
    losses = -c*l_0.*ones(N,1);
    U = wave_solver.make_step(-1i*c*n32,-1i*c*n32_t,losses,dt);    
    %set the boundaries... and obtain the final solution
    U = wave_solver.set_bdry(U(end),'no');
    %%%%%%%%%%%%%%%%%%%%%%%%   
    r11 = r11_solver.get_latest_solution();
    r33 = r33_solver.get_latest_solution();
    r22 = r22_solver.get_latest_solution();
    
    for p = 1:length(global_idx_rest)
        populations(:,p) = rate_eqn_solvers{p}.get_latest_solution();
    end
    
    r13 = r13_solver.get_latest_solution();
    n32 = n32_solver.get_latest_solution();
    n12 = n12_solver.get_latest_solution();
    
    t = t+dt;
    iter_ctr = iter_ctr + 1;
end
savename = [settings.name '_' settings.scenario '_NO_DISP_COMP_N_' num2str(settings.N) '_RING_RT_' num2str(settings.simRT)];
save(savename);

%get the current simulation parameters settings  
function [savename] = runsim_FP_RNFD( settingsfile,varargin )


close all;
lvar = length(varargin);


%handle user input!
if(length(varargin)>0)
    
    if mod(lvar,2) ~=0
        display('incorrect number of input arguments. Aborting!');
        return;
    end
    
    for i = 1:2:lvar
        name = varargin{i};
        val  = varargin{i+1};
        name = strtrim(lower(name)); %change to lower case!
        if strcmp(name,'init')
            init = val;
        else if strcmp(name,'scenario')
                scenario = (val)
            else if strcmp(name,'workspace')
                    %regexp to load all variables except init!
                    load(val,'-regexp','^(?!.*init.*).*$');
                else
                    display( ['Unknown input option name: ' name '. Aborting!' ]);
                    return;
                end
            end
        end
    end
else
    init = 1;
    scenario =1;
end

%parse all input settings params
settings = parseInput_fewlvl(settingsfile);

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

%%%% energies: E1 > E3 > E2 Optical transition is 3->2 TUNNELING transition is 1->3
E1 = settings.E1/hbar;      % in rad/ps
E3 = settings.E3/hbar;      % in rad/ps
E2 = settings.E2/hbar;      % in rad/ps
O13 = settings.O13/hbar;    % in rad/ps

E13 = E1-E3; %rad/ps; 1->3 traisition freq
E12 = E1-E2; %rad/ps; 1->2 transition freq
E32 = E3-E2; %rad/ps  3->2 transition freq (optical transition)
% central frequency and wave number!
E0 =  (E32+E12)/2; % central OPTICAL frequency.

%gain recovery times
INJ = settings.INJ; ULL = settings.ULL; LLL = settings.LLL; DEPOP = settings.DEPOP;

W = zeros(4);
for i = 1:4
    W(INJ,i) = settings.W_inj(i); W(ULL,i) = settings.W_ull(i);
    W(LLL,i) = settings.W_lll(i); W(DEPOP,i) = settings.W_dep(i);   
end

%%%% NOTICE THE CHANGE IN RATES HERE <- we do not include the backscattering rates w13 and w12!!! 
W31 = W(ULL,DEPOP) + W(ULL,INJ); W13 = W(INJ,ULL);
W21 = W(LLL,DEPOP) + W(LLL,INJ); W12 = W(INJ,LLL); 
W23 = W(LLL,ULL);                W32 = W(ULL,LLL);


G1 = W13 + W12;  G3 = W31 + W32; G2 = W21 + W23;

% % % % Added Pure Dephasing 
if(settings.deph>0)
    
    gamma_13 = 0.5*(G1+G3)+1/settings.Tdeph; %% dephsing of the resonant tunneling transition 
    gamma_32 = 0.5*(G2+G3)+1/settings.Tdeph; % dephasing of the optical transision...  
    gamma_12 = 0.5*(G2+G1)+1/settings.Tdeph; % dephasing of the latest transition 

else    
    gamma_13 = 0.5*(G1+G3); %% dephsing of the resonant tunneling transition
    gamma_32 = 0.5*(G2+G3); % dephasing of the optical transision...  
    gamma_12 = 0.5*(G2+G1); % dephasing of the latest transition 
end

dE13 = -1i*E13 - gamma_13; %
dE32 = 1i*(E0 - E32) - gamma_32; %
dE12 = 1i*(E0 - E12)- gamma_12; %


f0 = E0/2/pi;

if (scenario ==2)
    shb = 1;
else
    shb = -1;
end



Ncarriers = settings.dN*(100^3)*settings.Ld/settings.Lp; % cm^-3 --> m^-3; carrier density
Overlap = settings.Overlap;  % overlap factor -> dimensionless

%calculate the normalization constant, i.e. the number on which we divide
%all quantities ( except the electric field envelope) to simplify
% our initial system. it is important as it determines the initial value of
% the overall electron population inside the system-> a quantity that shall
% be perserved throughout the whole simulaiton!!

trace_rho = ((E0*1E12*Ncarriers*Overlap*((zUL*1E-9*Constants('q0'))^2))/(Constants('eps0')*n*Constants('c')*Constants('hbar')))/(1/(lch*tch));
%cavity loss l_0 in (cm^-1) --> l_0*100 in (m^-1) --> 1 mm^-1
l_0 = settings.loss*100/(1/lch);
%%%%%%%%%%%%%%%%% simuation parameters %%%%%%%%%%%%%%%%%
%grid size in x direction
N = settings.N; %nr of points in x direction
x = linspace(0,Ltot,N)';
dx = x(2) - x(1); dt = dx/c;
k0 = E0/c; D = settings.D*10^2/(1/tch); diffusion = 4*k0^2*D;

%%%% specify some of the main loop control parameters %%%%
iter_per_rt = round(T_R/dt);

tEnd = settings.simRT*T_R; % end time in tps
plotCtr = settings.plotCtr; %% set on how many iterations should the program plot the intensity
pulse_len = settings.recordRT*T_R; % how many ps should we record the pulse for

if(init > 0)
    clc; 
    t = 0;
    
    U = zeros(N,1); 
    V = 0*U;
    
    % population vectors together with first derivative vectors
    % Make sure that sum r_ii = trace_rho!
    r110 = trace_rho*(ones(N,1)*1); r110(1) = r110(end); 
    r330 = trace_rho*(ones(N,1)*0); r330(1) = r330(end); 
    r220 = (trace_rho-r110-r330); r220(1)=r220(end); 
    r11p = zeros(N,1);r33p = zeros(N,1);r22p = zeros(N,1);
    
    % coherence terms together with first derivative vectors!
    r130 =1E-15*((ones(N,1)-0.5) +1i*(ones(N,1)-0.5)); r130(1) = r130(end); 
    r13p = zeros(N,1); r13m = zeros(N,1);
    
    n32p = 1E-15*((rand(N,1)-0.5) +1i*(rand(N,1)-0.5)); n32p(1) = n32p(end); 
    n32m = 1E-15*((rand(N,1)-0.5) +1i*(rand(N,1)-0.5)); n32m(1) = n32m(end); 
    n12p = 1E-15*((rand(N,1)-0.5) +1i*(rand(N,1)-0.5)); n12p(1) = n12p(end); 
    n12m = 1E-15*((rand(N,1)-0.5) +1i*(rand(N,1)-0.5)); n12m(1) = n12m(end); 
    
    t = dt; nr_steps = settings.nr_steps;    
    % variable with changing lenght that stores the instantaneouse intensity at a predefined point inside the active region.
    I_t = 0; E_p = 0; E_m = 0 ; 
    r11_t =0;  r33_t =0;  r22_t =0; 
    pop_sum = 0; %variable with changing length that stores the instantaneous total population at idx!
    idx = 1; % the index of the spatial point we store the time-resolved simulated quantities at!
    iter_ctr = 0;  ctr = 1;   recordIter = 1 ; % nr of iterations after which to record new statistics!

    %initialize the solvers!
    
    r110_solver = MS(nr_steps,N,[],r110); r11p_solver =MS(nr_steps,N,[],r11p); 
    r330_solver = MS(nr_steps,N,[],r330); r33p_solver = MS(nr_steps,N,[],r33p); 
    r220_solver = MS(nr_steps,N,[],r220); r22p_solver = MS(nr_steps,N,[],r22p); 
    
    n12p_solver = MS(nr_steps,N,[],n12p); n12m_solver = MS(nr_steps,N,[],n12m);
    n32p_solver = MS(nr_steps,N,[],n32p); n32m_solver = MS(nr_steps,N,[],n32m); 
    
    r130_solver = MS(nr_steps,N,[],r130);
    r13p_solver = MS(nr_steps,N,[],r13p); r13m_solver = MS(nr_steps,N,[],r13m); 
    
    U_solver = RNFDSolver(N,dx,+1,c, U);
    V_solver = RNFDSolver(N,dx,-1,c,V);
    
end

while(t< tEnd)

    %%plot some of the results if neeed ariseth :D
    if(mod(iter_ctr,plotCtr) == 0)
        clc;
        display([settings.name ' simulation!']); 
        display(['iteration Nr. = ' num2str(iter_ctr) ' @ RT = ' num2str(t/T_R)])
        intensity = U.*conj(U) + V.*conj(V); maxInt = max(intensity);
        display(['max Intensity: ' num2str(maxInt) ]);
        ax = plotyy(x,intensity,x,[r110,r330,r220]);
        set(ax(1),'Xlim',[0 Ltot]);       
        set(ax(2),'Xlim',[0 Ltot]);       
        set(ax(2),'Ylim',[0 trace_rho]);
        title([ settings.name ' (MS + RNFD)  @ t = ' num2str(t)]);
        getframe;
        
        %if the intensity has decreased below 1E-43 after 1000 iterations
        %then abort!
        if(maxInt  < 1E-43 && iter_ctr > 1E3)
            display('LASING OFF!! ABORTING!')
            break;
        end
    end
    
    
    %%%% obtain the field, field intensity and the total population at
    %%%% position "idx" ...
    if(t >= tEnd - pulse_len && mod(iter_ctr,recordIter) == 0)
        
        intensity = U(idx).*conj(U(idx)) + V(idx).*conj(V(idx));
        
        I_t(ctr) = intensity;
        
        %forward and backward fields!
        E_p(ctr)= U(idx);
        E_m(ctr)= V(idx); 
        
        %store the populations! 
        r11_t(ctr)= r110(idx);
        r33_t(ctr)= r330(idx);
        r22_t(ctr) = r220(idx);
        pop_sum(ctr) = r110(idx)+r220(idx)+r330(idx);
        ctr = ctr+1;   
        
    end
    
    %%%%%% BLOCH PART %%%%%%
    %%%% POPULATIONS
    r110_t = 1i*O13.*(r130-conj(r130)) + r330*W31 + r220*W21 - G1*r110;
    r110_solver.make_step(r110_t,dt);
    
    lmInteraction = conj(U).*n32p + conj(V).*n32m;
    r330_t = 1i*O13.*(conj(r130) - r130) +1i/2.*(lmInteraction-conj(lmInteraction)) + r110*W13 + r220*W23 - G3*r330;
    r330_solver.make_step(r330_t,dt);
    
    r220_t = -1i/2.*(lmInteraction-conj(lmInteraction)) + r110*W12 + r330*W32 - G2*r220;
    r220_solver.make_step(r220_t,dt);
   
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
    
    if(shb > 0 )
        r11p_t = 1i*O13.*(r13p-conj(r13m)) + r33p*W31 + r22p*W21 - (G1+diffusion)*r11p;
        r11p_solver.make_step(r11p_t,dt);
        
        r33p_t = 1i*O13.*(conj(r13m)-r13p)+1i/2*(conj(V).*(n32p) -(U).*conj(n32m)) +  r11p*W13 + r22p*W23 - (G3+diffusion)*r33p;
        r33p_solver.make_step(r33p_t,dt);
        
        r22p_t = -1i/2*(conj(V).*(n32p) -(U).*conj(n32m)) +  r11p*W12 + r33p*W32 - (G2+diffusion)*r22p;
        r22p_solver.make_step(r22p_t,dt);    
        
        r31p_t = dE13*r13p + 1i*O13*(r11p-r33p) +1i/2*(conj(V).*n12p);
        r13p_solver.make_step(r31p_t,dt);
        
        r31m_t = dE13*r13m + 1i*O13*conj(r11p-r33p) + 1i/2*(conj(U).*n12m);
        r13m_solver.make_step(r31m_t,dt);
   
    end
    
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
    
    r130 = r130_solver.get_latest_solution(); n32p = n32p_solver.get_latest_solution();
    n32m = n32m_solver.get_latest_solution(); n12p = n12p_solver.get_latest_solution(); n12m = n12m_solver.get_latest_solution();
    
    if(shb > 0)
    
        r11p = r11p_solver.get_latest_solution();
        r33p = r33p_solver.get_latest_solution();
        r22p = r22p_solver.get_latest_solution();
        
        r13p = r13p_solver.get_latest_solution(); 
        r13m = r13m_solver.get_latest_solution();
   
    end
    
    t = t+dt;
    iter_ctr = iter_ctr + 1;

end
savename = [settings.name '_' num2str(scenario)];
save(savename);
end

folders = {'TB\MIT'};
files = dir(folders{1});

for n = 1:length(files)
    if( (strfind(files(n).name,'Input03')))
            inputfile = [folders{1} '\' files(n).name];
    end        
end

settingsfile = inputfile; 
fitting_file = 'fit_fin';
arguments = {'init',1,'scenario',2};


load(fitting_file);

close all;
lvar = length(arguments);


%handle user input! 
if(length(arguments)>0)
    
    if mod(lvar,2) ~=0
        display('incorrect number of input arguments. Aborting!');
        return;
    end
    
    for i = 1:2:lvar
        name = arguments{i};
        val  = arguments{i+1};
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

%%%
%%% CUSTOM MODIFIED : OFF
%%%

%parse all input files
settings = parseInputvII(settingsfile);

tch = settings.tch; % characteristic time..
lch = settings.lch; % characteristic length...

%speed of light in vacuum: 299792458m/s --> in tch/lch
c_0 = Constants('c',{'time',tch},{'length',lch});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEVICE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%
%cavity length  -->
Ltot = settings.Ltot; % mm
%phase velocity inside the medium ( in mm per picosecond ... )
n_eff = settings.n;  c = c_0/n_eff; T_R = 2*Ltot/c; f_R = 1/T_R;
%%%%dipole mtx elements (in C nm)
zUL = settings.zUL;
% hbar in eV-ps
hbar = Constants('hbar',{'time',tch})/Constants('q0');

% these four indices specify the core levels of our simulation, and namely
% the injector the ull , the lll and the depopulation level. All other
% possible levels will be included as rate equations in our simulation
INJ = settings.INJ; ULL = settings.ULL; LLL = settings.LLL; DEPOP = settings.DEPOP;


% % central frequency and wave number!
E0 =  3.8*2*pi;%(E32+E12)/2; % central OPTICAL frequency.



%extract the indices of the rest of the laser levels, which will be
%included in the system as simple rate equations

NLVLS = 7; % nr of levels to consider !
idx_rest = setdiff(1:NLVLS,[INJ,ULL,LLL,DEPOP]);
N_rest = length(idx_rest); %take the number of residual levels!
G = zeros(settings.N,NLVLS);

if (scenario ==2 || scenario == 4)
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
% be perserved throughout the whole simulaiton ! !

trace_rho = ((E0*1E12*Ncarriers*Overlap*((zUL*1E-9*Constants('q0'))^2))/(Constants('eps0')*n_eff*Constants('c')*Constants('hbar')))/(1/(lch*tch));


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
    
    x0 = Ltot/3;  tp = 6.0; aE_in = @(z,time) exp(-(time-(z-x0)/c).^2/tp^2);
    f_p = aE_in(x,0); f_p = 3*f_p/max(abs(f_p));  f_m = 0*f_p;
    
    %debye polarizations
    P_u = 0*f_p; P_v = 0*f_m; dP_u = P_u; dP_v = P_v;
    
    
    % population vectors together with first derivative vectors
    r110 = trace_rho*(ones(N,1)*0); r110(1) = r110(end);
    r330 = trace_rho*(ones(N,1)*1); r330(1) = r330(end);
    r220 = trace_rho*(ones(N,1)*0); r220(1)=r220(end);
    
    r11p = zeros(N,1);r33p = zeros(N,1);r22p = zeros(N,1);
    %rate populations
    populations = zeros(N,N_rest);
    
    % coherence terms together with first derivative vectors!
    r310 =0*((ones(N,1)-0.5) +1i*(ones(N,1)-0.5)); r310(1) = r310(end);
    r31p = zeros(N,1); r31m = zeros(N,1);
    
    n23p =  1E-15*((rand(N,1)-0.5) +1i*(rand(N,1)-0.5)); n23p(1) = n23p(end);
    n23m =  1E-15*((rand(N,1)-0.5) +1i*(rand(N,1)-0.5)); n23m(1) = n23m(end);
    
    n21p =  1E-15*((rand(N,1)-0.5) +1i*(rand(N,1)-0.5)); n21p(1) = n21p(end);
    n21m =  1E-15*((rand(N,1)-0.5) +1i*(rand(N,1)-0.5)); n21m(1) = n21m(end);
    
    t = dt; nr_steps = settings.nr_steps;
    
    I_t = 0; % variable with changing lenght that stores the instantaneouse intensity at a predefined point inside the active region.
    E_p = 0 ; E_m = 0 ;
    V_TL_t = 0; J_TL_t = 0; I_TL_t = 0; 
    r11_t =0;  r33_t =0;  r22_t =0;
    
    pop_sum = 0; %variable with changing length that stores the instantaneous total population at idx!
    idx = 1; % the index of the predefined point we refered to !
    iter_ctr = 0;    ctr = 1;   recordIter = 1 ; % nr of iterations after which to record new statistics!
    
    r110_solver = MS(nr_steps,N,[],r110); r11p_solver =MS(nr_steps,N,[],r11p);
    r330_solver = MS(nr_steps,N,[],r330); r33p_solver = MS(nr_steps,N,[],r33p);
    r220_solver = MS(nr_steps,N,[],r220); r22p_solver = MS(nr_steps,N,[],r22p);
    
    
    rate_eqn_solvers = [];
    for p=1:N_rest
        rate_eqn_solvers{p} =  MS(nr_steps,N,[],populations(:,p));
    end
    
    pu_solver =  MS(nr_steps,N,[],P_u); pv_solver =  MS(nr_steps,N,[],P_v);
    
    n21p_solver = MS(nr_steps,N,[],n21p); n21m_solver = MS(nr_steps,N,[],n21m);
    n23p_solver = MS(nr_steps,N,[],n23p); n23m_solver = MS(nr_steps,N,[],n23m);
    
    r310_solver = MS(nr_steps,N,[],r310);
    r31p_solver = MS(nr_steps,N,[],r31p); r31m_solver = MS(nr_steps,N,[],r31m);
    
    U_solver = RNFDSolver(N,dx,+1,c, f_p);
    V_solver = RNFDSolver(N,dx,-1,c,f_m);
    
    % transmission line params:
    bias = 10.2/1E1; %initial bias in kV/mm; 
    v_0 = bias;
    v_TL = v_0*ones(N,1); % transmission line voltage per unit length (in units kV/mm);
    J_TL = zeros(N,1); % QCL curr density (in units A/mm^2)
  
    %%% curr density calculation:
    rates = (r110).*W_fit{INJ,DEPOP}(v_TL) + (r220).*W_fit{LLL,DEPOP}(v_TL) + (r330).*W_fit{ULL,DEPOP}(v_TL);
    for p = 1:N_rest
        p_glob_idx = idx_rest(p);
        rates = rates + (populations(:,p)).*W_fit{p_glob_idx,DEPOP}(v_TL);
    end
    
    %convert the curr density in A/m^2 into A/mm^2 by dividing by 1E6;
    J_TL = (Constants('q0')*(settings.Lp*1E-9)*Ncarriers/trace_rho*1E12*rates)/1E6; %in A/mm^2
    eps_ch = 1E15*Constants('eps0');    mu_ch  = 1E3*Constants('mu0'); 
   %impedance:
    Z_0 = sqrt(mu_ch/eps_ch)/n_eff; 
    i_0 = Ltot*12; %trapz(x,J_TL); % v_0/Z_0;
    i_TL = i_0*ones(N,1);%(Ltot-x)./Ltot;

end

legendinfo = {'Intensity','inj','ull','lll'}; nlegend = length(legendinfo);
for p = 1:length(idx_rest)
    legendinfo{nlegend+p} = [' lvl: ' num2str(idx_rest(p)) ];
end

while(t< tEnd)
    
    
    %%%%%%%% transmission line equations %%%%%%%%%%
%     
%     %%% curr density calculation:
    rates = (r110).*W_fit{INJ,DEPOP}(v_TL) + (r220).*W_fit{LLL,DEPOP}(v_TL) + (r330).*W_fit{ULL,DEPOP}(v_TL);
    for p = 1:N_rest
        p_glob_idx = idx_rest(p);
        rates = rates + (populations(:,p)).*W_fit{p_glob_idx,DEPOP}(v_TL);
    end
    %convert the curr density in A/m^2 into A/mm^2 by dividing by 1E6;
    J_TL = (Constants('q0')*(settings.Lp*1E-9)*Ncarriers/trace_rho*1E12*rates)/1E6; %in A/mm^2
    
    i_TL(2:end) = i_TL(2:end)-(v_TL(2:end)-v_TL(1:end-1))/Z_0; i_TL(1) = i_0; 
    v_TL(1:end-1) = v_TL(1:end-1)-dt/(eps_ch*n_eff^2)*J_TL(1:end-1)-Z_0*(i_TL(2:end)-i_TL(1:end-1));
 
    %set bdry conditions for v at right end! 
    v_TL(end) = v_TL(end-1);
    %% put a constraint !!! ARTIFICIAL !!! on v_TL to be in the region 9-12 kV/cm (0.9-1.2 kv/mm)
    
    v_TL = (v_TL >= 0.9).*(v_TL<=1.2).*v_TL+ (v_TL < 0.9)*0.9 + (v_TL > 1.2)*1.2;
    
    %%%%% correctly setup the energies, the anticrossings and the scattering rates... 
    
    %obtain the energies of the core levels for the simulation
%     E1 = E_fit{INJ}(v_TL)/hbar; %resonance
    E3 = E_fit{ULL}(v_TL)/hbar; 
    E1 = E3; 
    E2 = E_fit{LLL}(v_TL)/hbar; % in rad/ps
    %%%% energies: E1 > E3 > E2 Optical transition is 3->2 TUNNELING transition
    %%%% is 1->3
    O13 = O_13_fit(v_TL)/hbar; % in rad/ps
    E13 = E1-E3; %rad/ps; 1->3 traisition freq
    E12 = E1-E2; %rad/ps; 1->2 transition freq
    E32 = E3-E2; %rad/ps  3->2 transition freq (optical transition)

    %extract the indices of the rest of the laser levels, which will be
    %included in the system as simple rate equations
    idx_rest = setdiff(1:NLVLS,[INJ,ULL,LLL,DEPOP]);
    N_rest = length(idx_rest); %take the number of residual levels!
    G = zeros(settings.N,NLVLS);
    
    for el = setdiff(1:NLVLS,[INJ,DEPOP])
        G(:,INJ) = G(:,INJ)+W_fit{INJ,el}(v_TL);
    end
  
    for lvl = setdiff(1:NLVLS,[INJ,DEPOP])
      for el = 1:NLVLS
        G(:,lvl) = G(:,lvl) +  W_fit{lvl,el}(v_TL);
      end
    end

    % % % % Added Pure Dephasing
    if(settings.deph>0)

        gamma_31 = 0.5*(G(:,INJ)+G(:,ULL)); %% dephsing of the resonant tunneling transition
        gamma_23 = 0.5*(G(:,LLL)+G(:,ULL))+1/settings.Tdeph; % dephasing of the optical transision...
        gamma_21 = 0.5*(G(:,LLL)+G(:,INJ))+1/settings.Tdeph; % dephasing of the latest transition

    else
        gamma_31 = 0.5*(G(:,INJ)+G(:,ULL)); %% dephsing of the resonant tunneling transition
        gamma_23 = 0.5*(G(:,LLL)+G(:,ULL)); % dephasing of the optical transision...
        gamma_21 = 0.5*(G(:,LLL)+G(:,INJ)); % dephasing of the latest transition
    end

    dE31 = 1i*E13 - gamma_31; %
    dE23 = -1i*(E0 - E32) - gamma_23; %
    dE21 = -1i*(E0 - E12)- gamma_21; %
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%plot some of the results if neeed ariseth :D
    
    if(mod(iter_ctr,plotCtr) == 0)
        clc;
        display([settings.name ' simulation!']);
        display(['iteration Nr. = ' num2str(iter_ctr) ' @ RT = ' num2str(t/T_R)])
        intensity = f_p.*conj(f_p) + f_m.*conj(f_m); maxInt = max(intensity);
        display(['max Intensity: ' num2str(maxInt) ]);
       
        subplot(3,1,1)
        ax = plotyy(x,intensity,x,[r110,r330,r220,populations]);
        title([' t = ' num2str(t)]);
        
        legend(legendinfo);
        set(ax(1),'Xlim',[0 Ltot]);
        set(ax(2),'Xlim',[0 Ltot]);
        set(ax(2),'Ylim',[0 trace_rho]);
        
        subplot(3,1,2) ;
        plotyy(x,v_TL*10,x,i_TL); legend('Bias (kv/cm)','TL current (A/mm)');
        %plot qcl curr density 
        subplot(3,1,3);
        plot(x,J_TL); xlabel('x direction [mm]'); ylabel('QCL curr density A/mm^2');
        
        getframe;
        
        if(maxInt  < 1E-43 && iter_ctr > 1E3)
            display('LASING OFF!! ABORTING!')
            break;
        end
    end
    
    
    %%%% obtain the field, field intensity and the total population at
    %%%% position "idx" ...
    if(t >= tEnd - pulse_len && mod(iter_ctr,recordIter) == 0)
        
        intensity = f_p(idx).*conj(f_p(idx)) + f_m(idx).*conj(f_m(idx));
        I_t(ctr) = intensity;
        E_p(ctr)= f_p(idx);
        E_m(ctr)= f_m(idx);
        V_TL_t(ctr) = v_TL(idx);
        I_TL_t(ctr) = i_TL(idx);
        J_TL_t(ctr) = J_TL(idx);
        

        r11_t(ctr)= r110(idx);
        r33_t(ctr)= r330(idx);
        r22_t(ctr) = r220(idx);
        pop_sum(ctr) = r110(idx)+r220(idx)+r330(idx);
        
        for p = 1:N_rest
            pop_sum(ctr) = pop_sum(ctr)+populations(idx,p);
        end
        ctr = ctr+1;
        
    end
    
    
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%% BLOCH PART %%%%%%
    %%%% POPULATIONS
    r110_t = 1i*O13.*(conj(r310)-r310) +(W_fit{ULL,INJ}(v_TL) + W_fit{ULL,DEPOP}(v_TL)).*r330 + ( W_fit{LLL,INJ}(v_TL)+W_fit{LLL,DEPOP}(v_TL) ).*r220 - G(:,INJ).*r110;
    
    for p = 1:N_rest
        p_glob_idx = idx_rest(p);
        r110_t = r110_t + (W_fit{p_glob_idx,INJ}(v_TL)+W_fit{p_glob_idx,DEPOP}(v_TL)).*populations(:,p);
    end
    r110_solver.make_step(r110_t,dt);
    
    lmInteraction = f_p.*conj(n23p) + f_m.*conj(n23m);
    r330_t = 1i*O13.*(r310 - conj(r310)) +1i/2.*(lmInteraction-conj(lmInteraction)) + r110.*W_fit{INJ,ULL}(v_TL) + r220.*W_fit{LLL,ULL}(v_TL) - G(:,ULL).*r330;
    for p = 1:N_rest
        p_glob_idx = idx_rest(p);
        r330_t = r330_t + W_fit{p_glob_idx,ULL}(v_TL).*populations(:,p);
    end
    r330_solver.make_step(r330_t,dt);
    
    r220_t = -1i/2.*(lmInteraction-conj(lmInteraction)) + r110.*W_fit{INJ,LLL}(v_TL) + r330.*W_fit{ULL,LLL}(v_TL) - G(:,LLL).*r220;
    for p = 1:N_rest
        p_glob_idx = idx_rest(p);
        r220_t = r220_t + W_fit{p_glob_idx,LLL}(v_TL).*populations(:,p);
    end
    r220_solver.make_step(r220_t,dt);
    
    for p = 1:N_rest
        
        p_glob_idx =idx_rest(p);
        rpp_0_t = W_fit{INJ,p_glob_idx}(v_TL).*r110+W_fit{ULL,p_glob_idx}(v_TL).*r330+W_fit{LLL,p_glob_idx}(v_TL).*r220 - G(:,p_glob_idx).*populations(:,p);
        
        %add the rate equations part!
        for j =1:N_rest % nr
            if j~= p
                j_glob_idx =idx_rest(j);
                rpp_0_t = rpp_0_t + W_fit{j_glob_idx,p_glob_idx}(v_TL).*populations(:,j);
            end
        end
        rate_eqn_solvers{p}.make_step(rpp_0_t,dt);
    end
    
    
    %%%% COHERENCES
    r310_t = dE31.*r310 - 1i*O13.*(r110-r330) -1i/2*(conj(f_p).*n21p + conj(f_m).*n21m);
    r310_solver.make_step(r310_t,dt);
    
    n23p_t = dE23.*n23p - 1i/2*(f_p.*(r330-r220) + f_m.*conj(r33p-r22p)) + 1i*O13.*n21p;
    n23p_solver.make_step(n23p_t,dt);
    
    n23m_t = dE23.*n23m - 1i/2*(f_m.*(r330-r220) + f_p.*(r33p-r22p)) + 1i*O13.*n21m;
    n23m_solver.make_step(n23m_t,dt);
    
    n21p_t = dE21.*n21p -1i/2*(f_p.*r310 + f_m.*r31m) + 1i*O13.*n23p;
    n21p_solver.make_step(n21p_t,dt);
    
    n21m_t = dE21.*n21m - 1i/2*(f_m.*r310 +f_p.*r31p) + 1i*O13.*n23m;
    n21m_solver.make_step(n21m_t,dt);
    
    if(shb > 0 )
        r11p_t = 1i*O13.*(conj(r31m)-r31p) + (W_fit{ULL,INJ}(v_TL)+W_fit{ULL,DEPOP}(v_TL)).*r33p+ (W_fit{LLL,INJ}(v_TL)+W_fit{LLL,DEPOP}(v_TL)).*r22p - (G(:,INJ)+diffusion).*r11p;
        r11p_solver.make_step(r11p_t,dt);
        
        r33p_t = 1i*O13.*(r31p-conj(r31m))+1i/2*(f_m.*conj(n23p) -conj(f_p).*n23m) +  W_fit{INJ,ULL}(v_TL).*r11p + W_fit{LLL,ULL}(v_TL).*r22p - (G(:,ULL)+diffusion).*r33p;
        r33p_solver.make_step(r33p_t,dt);
        
        r22p_t = -1i/2*(f_m.*conj(n23p) -conj(f_p).*n23m) +  W_fit{INJ,LLL}(v_TL).*r11p + W_fit{ULL,LLL}(v_TL).*r33p - (G(:,LLL)+diffusion).*r22p;
        r22p_solver.make_step(r22p_t,dt);
        
        r31p_t = dE31.*r31p - 1i*O13.*(r11p-r33p) -1i/2*(conj(f_p).*n21m);
        r31p_solver.make_step(r31p_t,dt);
        
        r31m_t = dE31.*r31m  -1i*O13.*conj(r11p-r33p) - 1i/2*(conj(f_m).*n21p);
        r31m_solver.make_step(r31m_t,dt);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%% Optical Field %%%%%%%%%%%%%%%%%%%%%%%%
    losses = -c*l_0.*ones(N,1);
    e_1 = 0;
    f_p = U_solver.make_step(1i*c*n23p+1i*e_1*c*P_u,1i*c*n23p_t+1i*e_1*c*dP_u,losses,dt);
    f_m = V_solver.make_step(1i*c*n23m+1i*e_1*c*P_v,1i*c*n23m_t+1i*e_1*c*dP_v,losses,dt);
    
    %set the boundaries... and obtain the final solution
    f_p = U_solver.set_bdry(f_m(1),'no');
    f_m = V_solver.set_bdry('no',f_p(N));
    %%%%%%%%%%%%%%%%%%%%%%%% ************ %%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    P_u = pu_solver.get_latest_solution();   P_v = pv_solver.get_latest_solution();
    r110 = r110_solver.get_latest_solution();
    r330 = r330_solver.get_latest_solution();
    r220 = r220_solver.get_latest_solution();
    
    for p = 1:length(idx_rest)
        populations(:,p) = rate_eqn_solvers{p}.get_latest_solution();
    end
    
    r310 = r310_solver.get_latest_solution(); n23p = n23p_solver.get_latest_solution();
    n23m = n23m_solver.get_latest_solution(); n21p = n21p_solver.get_latest_solution(); n21m = n21m_solver.get_latest_solution();
    if(shb > 0)
        r11p = r11p_solver.get_latest_solution();
        r33p = r33p_solver.get_latest_solution();
        r22p = r22p_solver.get_latest_solution();
        r31p = r31p_solver.get_latest_solution();
        r31m = r31m_solver.get_latest_solution();
    end
    t = t+dt;
    iter_ctr = iter_ctr + 1;
end
savename = [settings.name '_TL_' num2str(settings.simRT)];
save(savename);


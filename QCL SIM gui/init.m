function init(settings)

tch = settings.tch; % characteristic time.. seconds per picosecond
lch = settings.lch; % characteristic length... meters per millimeter

%speed of light in vacuum: 299792458m/s --> in tch/lch
c_0 = Constants('c',{'time',tch},{'length',lch});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEVICE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%
%cavity length  -->
Ltot = settings.; % mm

%phase velocity inside the medium ( in mm per picosecond ... )
f_R = 6.8E-3; %6.8 GHz
T_R = 1/f_R;  %roundtrip time
n = (T_R*c_0)/(2*Ltot);
c = c_0/n;

%%%%dipole mtx elements (in Cnm)
z23 = 3.00;
%%%%dipole coupling strengths ratio:
hbar = 6.58211928E-16/tch; % eV�ps

%%%% energies: E1 > E3 > E2 Optical transition is 3->2 TUNNELING transition
%%%% is 1->3
E1 = 215.13E-3/hbar;
E3 = 215.13E-3/hbar;
E2 = 199.43E-3/hbar; % in rad/ps
O13 = -1.0E-3/hbar; % in rad/ps

E13 = E1-E3; %rad/ps; 1->3 traisition freq
E12 = E1-E2; %rad/ps; 1->2 transition freq
E32 = E3-E2; %rad/ps  3->2 transition freq (optical transition)
% central frequency and wave number!
E0 = (E12 + E32)/2; % central OPTICAL frequency.


%gain recovery times
%%% FROM IVAN' MC calculation
w1p3 = 0.0178;   
w1p2 = 0.0024;
w1p1 = 0.0020;

w31p = 0.0126;
w32 = 0.0852;
w31 = 0.0497;

w21p = 0.0055;	
w23  = 0.0497;
w21  = 0.8087;

w11p = 0.0019;
% w13  = 0.0042;	
% w12  = 0.0037;	

W31 = w31 + w31p; W13 = w1p3; %%%% NOTICE THE CHANGE IN RATES HERE!!! 
W21 = w21 + w21p; W12 = w1p2; %%%% NOTICE THE CHANGE IN RATES HERE <- we do not include the backscattering rates w13 and w12!!! 
W23 = w23; W32 = w32;

G1 = W13 + W12; 
G3 = W31 + W32;
G2 = W21 + W23;

% Tdeph = 0.75; 
% % % % Added Pure Dephasing 
% gamma_31 = 0.5*(G1+G3)+1/Tdeph; %% dephsing of the resonant tunneling transition 
% gamma_23 = 0.5*(G2+G3)+1/Tdeph; % dephasing of the optical transision...  
% gamma_21 = 0.5*(G2+G1)+1/Tdeph; % dephasing of the latest transition 

% Tdeph = 0.75; 
% No Pure Dephasing!!!
gamma_31 = 0.5*(G1+G3); %% dephsing of the resonant tunneling transition 
gamma_23 = 0.5*(G2+G3); % dephasing of the optical transision...  
gamma_21 = 0.5*(G2+G1); % dephasing of the latest transition 

%vacuum permitivity
eps0 = 8.854187817E-12;% F/m
e0 = 1.60217657E-19; % C elementary charge
Ncarriers = 1E16*(100^3); % cm^-3 --> m^-3; carrier density
hbar_full = 1.054571726E-34; %J·s planck's reduced constant in Joules per sec.
Overlap = 0.8;  % overlap factor -> dimensionless

dE31 = 1i*E13 - gamma_31; %
dE23 = -1i*(E0 - E32) - gamma_23; %
dE21 = -1i*(E0 - E12)-gamma_21; %


k0 = (E0/c)/lch; % central mode's wavenumber in units 1/m
%calculate the normalization constant, i.e. the number on which we divide
%all quantities ( except the electric field envelope) to simplify
% our initial system. it is important as it determines the initial value of
% the overall electron population inside the system-> a quantity that shall
% be perserved throughout the whole simulaiton ! !

trace_rho = ((k0*Ncarriers*Overlap*((z23*1E-9*e0)^2))/(eps0*n^2*hbar_full))/(1/(lch*tch));

%cavity loss l = 10 (cm^-1) --> 1000m^-1 --> 1 mm^-1
l_0 = 1000/(1/lch);

%%%%%%%%%%%%%%%%% simuation parameters %%%%%%%%%%%%%%%%%


%grid size in x direction
N = 3000; %nr of points in x direction1
x = linspace(0,Ltot,N)';
dx = x(2) - x(1);
dt = dx/c;


k0_mm = E0/c;
D = 46*10^2/(1/tch);
diffusion = 4*k0_mm^2*D;

T_g = 1/(diffusion + G3);

%%%% specify some of the mainloop control parameters %%%%
iter_per_rt = round(T_R/dt);

tEnd = 1000*T_R; % end time in tps
plotCtr = 1000; %% set on how many iterations should the program plot the intensity
pulse_len = tEnd; % how many ps should we record the pulse for

if(init > 0)
    
    clc;
    t = 0;
    pField = zeros(N,1);
    mField = pField;
    % population vectors together with first derivative vectors
    r110 = trace_rho*(ones(N,1)*1); r110(1) = r110(end); 
    r330 = trace_rho*(ones(N,1)*0); r330(1) = r330(end); 
    r220 = (trace_rho-r110-r330); r220(1)=r220(end); 
    r11p = zeros(N,1);r33p = zeros(N,1);r22p = zeros(N,1);
    
    % coherence terms together with first derivative vectors!
    r310 =0*((ones(N,1)-0.5) +1i*(ones(N,1)-0.5)); r310(1) = r310(end); 
    r31p = zeros(N,1); r31m = zeros(N,1);
    
    n23p = 1E-6*((rand(N,1)-0.5) +1i*(rand(N,1)-0.5)); n23p(1) = n23p(end); 
    n23m = 1E-6*((rand(N,1)-0.5) +1i*(rand(N,1)-0.5)); n23m(1) = n23m(end); 
    
    n21p = 1E-6*((rand(N,1)-0.5) +1i*(rand(N,1)-0.5)); n21p(1) = n21p(end); 
    n21m = 1E-6*((rand(N,1)-0.5) +1i*(rand(N,1)-0.5)); n21m(1) = n21m(end); 
    
    t = dt;  
    nr_steps = 5;    
    
    l = l_0;
    I_t = 0; % variable with changing lenght that stores the instantaneouse intensity at a predefined point inside the active region.
    E_p = 0 ; 
    E_m = 0 ; 
    
    r11_t =0;
    r33_t =0; 
    r22_t =0; 
    
    pop_sum = 0; % variable with changing length that stores the instantaneous total population at idx!
    idx = 1; % the index of the predefined point we refered to !
    iter_ctr = 0;
    ctr = 1;
    recordIter = 1 ; % nr of iterations after which to record new statistics!

    r110_solver = MS(nr_steps,N,[],r110); r11p_solver =MS(nr_steps,N,[],r11p); 
    r330_solver = MS(nr_steps,N,[],r330); r33p_solver = MS(nr_steps,N,[],r33p); 
    r220_solver = MS(nr_steps,N,[],r220); r22p_solver = MS(nr_steps,N,[],r22p); 
    
    
    n21p_solver = MS(nr_steps,N,[],n21p); n21m_solver = MS(nr_steps,N,[],n21m);
    n23p_solver = MS(nr_steps,N,[],n23p); n23m_solver = MS(nr_steps,N,[],n23m); 
    
    r310_solver = MS(nr_steps,N,[],r310);
    r31p_solver = MS(nr_steps,N,[],r31p); r31m_solver = MS(nr_steps,N,[],r31m); 
    
    forward_solver =RNFDSolver(N,dx,+1,c, pField);
    backward_solver =RNFDSolver(N,dx,-1,c, mField);
end
clear;clc;close all;
sec1file = 'sec1.set';

interpDataFile = 'fitted_data1.mat';
%parse all input files and load the scatterin rates file !
params_s1 = input_parser(sec1file);

load(interpDataFile);
clc

bias = (params_s1.bias);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZE DEVICE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%
%phase velocity inside the medium ( in mm per picosecond ... )RT time and
%frequency
tch = params_s1.tch;
lch = params_s1.lch;
c = Constants('c',{'time',tch},{'length',lch})/params_s1.nTHz;

params_s1.c = c;
Ltot = params_s1.length;


T_R = 2*Ltot/c; f_R = 1/T_R;
params_s1.INJ = INJ; params_s1.ULL = ULL; params_s1.LLL = LLL; params_s1.DEPOP = DEPOP;
params_s1.NLVLS = NLVLS; params_s1.N_rest = params_s1.NLVLS-4;

hbar = Constants('hbar',{'time',tch})/Constants('q0');
%central field frequency.
E0 =  (E_fit{ULL}(bias/10)-E_fit{LLL}(bias/10))/hbar;

params_s1.E0 = E0; 

N = params_s1.N_pts;
%grid size in x direction
x = linspace(0,Ltot,N)';
dx = x(2) - x(1); 
dt = dx/c;


params_s1.dx = dx; params_s2.dx = dx;

params_s1.zUL = zUL_fit(bias/10);

params_s1.zNORM = params_s1.zUL;

params_s1.Ncarriers_cm = params_s1.dN*params_s1.Ld/params_s1.Lp; %carrier density

params_s1.IDX = 1:params_s1.N_pts;
params_s1.x = x(params_s1.IDX); 

record_U= 1; record_V = 1;

record_r110 = 1; record_r330 = 1; record_r220 = 1;
record_rRES = cell(NLVLS-4,1);
record_popsum = 1; record_v_TL =1;
record_i_TL = 1; record_J_TL =1;


TL_bias = bias/10*ones(N,1);

dm_model_s1 = DM_MODEL_3_LVL_RWA_FP(params_s1);
dm_model_s1.interpolate(TL_bias(params_s1.IDX),W_fit,E_fit,AC_fit,zUL_fit);

dat.N = N; dat.c = c; dat.dx = dx;
dat.dt = dt;
dat = makeMaxwellVars(dat);

rlgc1.C = 1;    %unit: pF/mm
rlgc1.L = 1.6e2;  %unit: pH/mm
rlgc2.C = 1;    %unit: pF/mm
rlgc2.L = 1.6e2;  %unit: pH/mm



%%%% specify some of the mainloop control parameters %%%%
idx = 1; ctr = 1; iter_ctr = 1;
iter_per_rt = round(T_R/dat.dt); 
dat.simRT = 100; tEnd = dat.simRT*T_R; % end time in tps


%simulation info storage arrays -> preallocate
recordingduration = tEnd; % how many ps should we record the pulse for
iterperrecord = 1; recordingiter  = round(recordingduration/iterperrecord/dt);
padsize = double(recordingiter-length(record_U));

t = 0;
P = zeros(dat.N,1); P_t = zeros(dat.N,1); M = P; M_t = P_t; losses = P_t;
interpCtr = 1; %set how often to interpolate the energies, scattering rates, dipole elements etc.
checkptIter = 1040000;% 1039800; %100000
f_plot = 1000;
f_display = 100;


while( t< tEnd)
    
    dm_model_s1.propagate(dat.U,dat.V,dt);
    
    [P1,P1_t,M1,M1_t,losses1] = dm_model_s1.get_polarization_and_losses();
    P(params_s1.IDX) = P1; P_t(params_s1.IDX) = P1_t;
    M(params_s1.IDX) = M1; M_t(params_s1.IDX) = M1_t;
    losses(params_s1.IDX) = losses1;
    
 

    dat = stepWave(dat,P,P_t,M,M_t,losses);

    dm_model_s1.update_state();


    if(mod(iter_ctr,f_display) == 0)
        clc;
        intensity = abs(dat.U).^2 + abs(dat.V).^2;
        display(['max intensity ' num2str(max(intensity))])
        plot(x,intensity);
        title(['t=',num2str(t),' ps','        |E_z|^2','        Iteration= ',num2str(iter_ctr)]);
        getframe;        

    end    
    
    t = t+dt;
    iter_ctr = iter_ctr + 1;
    
end




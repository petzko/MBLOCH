clear;clc;close all;
sec1file = 'sec1.set';
sec2file = 'sec2.set';

interpDataFile = 'fitted_data1.mat';
%parse all input files and load the scatterin rates file !
params_s1 = input_parser(sec1file);
params_s2 = input_parser(sec2file);

load(interpDataFile);
clc

bias = 1/2*(params_s1.bias + params_s2.bias);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZE DEVICE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%
%phase velocity inside the medium ( in mm per picosecond ... )RT time and
%frequency
tch = params_s1.tch;
lch = params_s1.lch;
c = Constants('c',{'time',tch},{'length',lch})/params_s1.nTHz;

params_s1.c = c;
params_s2.c = c;
Ltot = params_s1.length+params_s2.length;

T_R = 2*Ltot/c; f_R = 1/T_R;
params_s1.INJ = INJ; params_s1.ULL = ULL; params_s1.LLL = LLL; params_s1.DEPOP = DEPOP;
params_s1.IGNORELEVEL = 1; params_s1.NLVLS = NLVLS;

params_s2.INJ = INJ; params_s2.ULL = ULL; params_s2.LLL = LLL; params_s2.DEPOP = DEPOP;
params_s2.NLVLS = NLVLS; params_s2.IGNORELEVEL = 1;


hbar = Constants('hbar',{'time',tch})/Constants('q0');
%central field frequency.
E0 =  (E_fit{ULL}(bias/10)-E_fit{LLL}(bias/10))/hbar;

params_s1.E0 = E0;
params_s2.E0 = E0;

N = params_s1.N_pts+params_s2.N_pts;
%grid size in x direction
x = linspace(0,Ltot,N)';
dx = x(2) - x(1);
dt = dx/c;


params_s1.dx = dx; params_s2.dx = dx;
params_s1.dt = dt; params_s2.dt = dt;

params_s1.zUL = zUL_fit(bias/10);
params_s2.zUL = params_s1.zUL;

params_s1.zNORM = params_s1.zUL;
params_s2.zNORM = params_s2.zUL;

params_s1.Ncarriers_cm = params_s1.dN*params_s1.Ld/params_s1.Lp; %carrier density
params_s2.Ncarriers_cm = params_s2.dN*params_s2.Ld/params_s2.Lp; %carrier density

params_s1.IDX = 1:params_s1.N_pts;
params_s2.IDX = params_s1.N_pts+1:params_s1.N_pts+params_s2.N_pts;
params_s1.x = x(params_s1.IDX); params_s2.x = x(params_s2.IDX);

record_U= 1; record_V = 1;

record_r110 = 1; record_r330 = 1; record_r220 = 1;
record_rRES = cell(NLVLS-4,1);
record_popsum = 1; record_v_TL =1;
record_i_TL = 1; record_J_TL =1;


TL_bias = bias/10*ones(N,1);

dm_model_s1 = DM_MODEL_PSEUDO_3_LVL_RWA_FP(params_s1);
dm_model_s1.interpolate(TL_bias(params_s1.IDX),W_fit,E_fit,AC_fit,zUL_fit);

dm_model_s2 = DM_MODEL_PSEUDO_3_LVL_RWA_FP(params_s2);
dm_model_s2.interpolate(TL_bias(params_s2.IDX),W_fit,E_fit,AC_fit,zUL_fit);

dat.N = N; dat.c = c; dat.dx = dx;
dat.dt = dt;
dat = makeMaxwellVars(dat);

rlgc1.C = 2;    %unit: pF/mm
rlgc1.L = 1.6e2;  %unit: pH/mm
rlgc2.C = 2;    %unit: pF/mm
rlgc2.L = 1.6e2;  %unit: pH/mm

TL_model_s1 = TL_solver(params_s1,rlgc1);
TL_model_s2 = TL_solver(params_s2,rlgc2);

VV_TL1 = TL_model_s1.v_TL;
VV_TL2 = TL_model_s2.v_TL;
i1old = params_s2.current;


%%%% specify some of the mainloop control parameters %%%%
idx = 1; ctr = 1; iter_ctr = 1;
iter_per_rt = round(T_R/dat.dt);
dat.simRT = 100; tEnd = dat.simRT*T_R; % end time in tps


%simulation info storage arrays -> preallocate
recordingduration = tEnd; % how many ps should we record the pulse for
iterperrecord = 1; recordingiter  = round(recordingduration/iterperrecord/dt);
padsize = double(recordingiter-length(record_U));

dat.t = 0;
P = zeros(dat.N,1); P_t = zeros(dat.N,1); M = P; M_t = P_t; losses = P_t;
interpCtr = 1; %set how often to interpolate the energies, scattering rates, dipole elements etc.
checkptIter = 1040000;% 1039800; %100000
f_plot = 1000;
f_display = 100;

r11 = zeros(N,1); r22 = zeros(N,1); r33 = zeros(N,1);


while( dat.t< tEnd)
    
    dm_model_s1.propagate(dat.U,dat.V,dt);
    dm_model_s2.propagate(dat.U,dat.V,dt);
    
    [P1,P1_t,M1,M1_t,losses1] = dm_model_s1.get_polarization_and_losses();
    P(params_s1.IDX) = P1; P_t(params_s1.IDX) = P1_t;
    M(params_s1.IDX) = M1; M_t(params_s1.IDX) = M1_t;
    losses(params_s1.IDX) = losses1;
    
    
    [P2,P2_t,M2,M2_t,losses2] = dm_model_s2.get_polarization_and_losses();
    P(params_s2.IDX) = P2; P_t(params_s2.IDX) = P2_t;
    M(params_s2.IDX) = M2; M_t(params_s2.IDX) = M2_t;
    losses(params_s2.IDX) = losses2;
    
    J_TL1 = dm_model_s1.get_current_density(params_s1);
    J_TL2 = dm_model_s2.get_current_density(params_s2);
    J_tot1 = trapz(x(1:params_s1.N_pts),J_TL1);
    J_tot2 = trapz(x(1:params_s2.N_pts),J_TL2);
    
    dat = stepWave(dat,P,P_t,M,M_t,losses);
    
    %   Transmission line equations
    TL_model_s1.propagate_1(J_TL1);
    TL_model_s2.propagate_2(J_tot1,J_TL2);
    
    
    dm_model_s1.update_state();
    dm_model_s2.update_state();
    
    MM1 = TL_model_s1.v_TL > 1.22; NN1 = TL_model_s1.v_TL < 0.88;
    VV_TL1 = TL_model_s1.v_TL;
    VV_TL1(MM1)=1.22; VV_TL1(NN1)=0.88;
    
    MM2 = TL_model_s2.v_TL > 1.22; NN2 = TL_model_s2.v_TL < 0.88;
    VV_TL2 = TL_model_s2.v_TL;
    VV_TL2(MM2)=1.22; VV_TL2(NN2)=0.88;
    
    dm_model_s1.interpolate(VV_TL1,W_fit,E_fit,AC_fit,zUL_fit);%TL_model_s1.v_TL
    dm_model_s2.interpolate(VV_TL2,W_fit,E_fit,AC_fit,zUL_fit);%TL_model_s2.v_TL
    
    if(mod(iter_ctr,f_display) == 0)
        clc;
        [trace10,trace12] = dm_model_s1.get_avg_trace();
        [trace20,trace22] = dm_model_s2.get_avg_trace();
        display(['trace section 1: ' num2str(trace10) , '; ' num2str(trace12)]);
        display(['trace section 2: ' num2str(trace20) , '; ' num2str(trace22)]);
        display(['Iteration: ' num2str(iter_ctr)]);
        display(['average bias (kV/cm) sec 1: ' num2str(mean(TL_model_s1.v_TL)*10)]);
        display(['average bias (kV/cm) sec 2: ' num2str(mean(TL_model_s2.v_TL)*10)]);
        display(['v_1TL(1) (V): ' num2str(TL_model_s1.v_TL(1)*10)]);
        display(['v_2TL(1) (V): ' num2str(TL_model_s2.v_TL(1)*10)]);
        
        display(['i_in_1 (A): ' num2str(TL_model_s1.i_TL(1)*params_s1.width)]);
        display(['i_out_1 (A): ' num2str(J_tot1*params_s1.width)]);
        
        display(['i_in_2 (A): ' num2str(TL_model_s2.i_TL(1)*params_s1.width)]);
        display(['i_out_2 (A): ' num2str(J_tot2*params_s1.width)]);
        
        subplot(4,1,1);
        %
        r33(dm_model_s1.IDX) = dm_model_s1.r330;
        r22(dm_model_s1.IDX) = dm_model_s1.r220;
        r11(dm_model_s1.IDX) = dm_model_s1.r110;
        
        r33(dm_model_s2.IDX) = dm_model_s2.r330;
        r22(dm_model_s2.IDX) = dm_model_s2.r220;
        r11(dm_model_s2.IDX) = dm_model_s2.r110;
        
        
        
        plotyy(x,[r33,r22,r11],x,[abs(dat.U).^2,abs(dat.V).^2]);
        legend('\rho_{u}','\rho_{l}','\rho_{d}','|U|^2','|V|^2');
        
        title(['t=',num2str(dat.t),' ps',' Iteration= ',num2str(iter_ctr)]);
        
        subplot(4,1,2);
        J_TL(1:params_s1.N_pts) = J_TL1; 
        J_TL(params_s1.N_pts+1:params_s1.N_pts+params_s2.N_pts) = J_TL2;
        plot(x,J_TL);
        title('Current density J [A/mm^2]');
        
        subplot(4,1,3);
        V_TL(1:params_s1.N_pts) = TL_model_s1.v_TL; V_TL(params_s1.N_pts+1:params_s1.N_pts+params_s2.N_pts) = TL_model_s2.v_TL;
        plot(x,V_TL*params_s1.height*1e3);
        title('V [V]');
        
        subplot(4,1,4);
        I_TL(1:params_s1.N_pts) = TL_model_s1.i_TL; I_TL(params_s1.N_pts+1:params_s1.N_pts+params_s2.N_pts) = TL_model_s2.i_TL;
        plot(x,I_TL*params_s1.width);
        title('I [A]');
        getframe;
    end
    
    
    
    
    %%%% obtain the field, field intensity and the total population at position "idx" ...
    %store fields info
    record_U(iter_ctr)= dat.U(idx);  record_V(ctr)= dat.V(idx);
    
    dat.t = dat.t+dt;
    iter_ctr = iter_ctr + 1;
    
end

mydft(record_U,dt);

simname = 'sim-test';

savename = [simname];
save(savename);


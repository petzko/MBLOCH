function acmodelocking(simfilename,interpDataFile, savename)

%parse all input files and load the scatterin rates file !
params_s1 = input_parser(simfilename);
load(interpDataFile);


bias = params_s1.bias ;


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
params_s1.IGNORELEVEL = -1; params_s1.NLVLS = NLVLS;



hbar = Constants('hbar',{'time',tch})/Constants('q0');
%central field frequency.
E0 =  (E_fit{ULL}(bias/10)-E_fit{LLL}(bias/10))/hbar;

params_s1.E0 = E0;

N = params_s1.N_pts;
%grid size in x direction
x = linspace(0,Ltot,N)';
dx = x(2) - x(1);
dt = dx/c;


params_s1.dx = dx;
params_s1.dt = dt;

params_s1.zUL = zUL_fit(bias/10);

params_s1.zNORM = params_s1.zUL;

params_s1.Ncarriers_cm = params_s1.dN*params_s1.Ld/params_s1.Lp; %carrier density

params_s1.IDX = 1:params_s1.N_pts;
params_s1.x = x(params_s1.IDX);

record_U= 1; 

TL_bias = bias/10*ones(N,1);

dm_model_s1 = DM_MODEL_PSEUDO_3_LVL_RWA_FP(params_s1);
dm_model_s1.interpolate(TL_bias(params_s1.IDX),W_fit,E_fit,AC_fit,zUL_fit);


dat.N = N; dat.c = c; dat.dx = dx;
dat.dt = dt;
dat = makeMaxwellVars(dat);

rlgc.C = 1;     %unit: pF/mm
rlgc.L = 1.6e2; %unit: pH/mm
rlgc.R = 2;     %unit: Ohm/mm

TL_model_s1 = TL_solver(params_s1,rlgc);



%%%% specify some of the mainloop control parameters %%%%
idx = 500; iter_ctr = 1;
dat.simRT = 200; tEnd = dat.simRT*T_R; % end time in tps

%simulation info storage arrays -> preallocate
recordingduration = tEnd; % how many ps should we record the pulse for
iterperrecord = 1; recordingiter  = round(recordingduration/iterperrecord/dt);
padsize = double(recordingiter-length(record_U));

record_U= zeros(padsize,1);
record_V = zeros(padsize,1);
record_v_TL =zeros(padsize,1);
record_i_TL = zeros(padsize,1); 

record_r110 = zeros(padsize,1); 
record_r330 = zeros(padsize,1); 
record_r220 = zeros(padsize,1);
record_rRES = cell(NLVLS-4,padsize);
record_popsum = zeros(padsize,1);
record_J_TL =zeros(padsize,1);


dat.t = 0;
P = zeros(dat.N,1); P_t = zeros(dat.N,1); M = P; M_t = P_t; losses = P_t;
f_display = 500;



while( dat.t< tEnd)

    dm_model_s1.propagate(dat.U,dat.V,dt);

    [P1,P1_t,M1,M1_t,losses1] = dm_model_s1.get_polarization_and_losses();
    P(params_s1.IDX) = P1; P_t(params_s1.IDX) = P1_t;
    M(params_s1.IDX) = M1; M_t(params_s1.IDX) = M1_t;
    losses(params_s1.IDX) = losses1;


    J_TL1 = dm_model_s1.get_current_density(params_s1);
    J_tot1 = trapz(x(1:params_s1.N_pts),J_TL1);


    dat = stepWave(dat,P,P_t,M,M_t,losses);

    %   Transmission line equations
    if iter_ctr >= 2000
        if iter_ctr == 2000        %Set initial current distribution
           for mm = 1: N-1
           TL_model_s1.i_TL(mm) = J_tot1*(N-mm)/N;
           end 
        end
       
        TL_model_s1.propagate(J_TL1,f_R,(dat.t-dt*2000))
    end


    dm_model_s1.update_state();


%     MM1 = TL_model_s1.v_TL > 1.4; NN1 = TL_model_s1.v_TL < 0.7;
%     VV_TL1 = TL_model_s1.v_TL;
%     VV_TL1(MM1)=1.4; VV_TL1(NN1)=0.7;

    dm_model_s1.interpolate(TL_model_s1.v_TL,W_fit,E_fit,AC_fit,zUL_fit);


    if(mod(iter_ctr,f_display) == 0)
        clc;
        [trace10,trace12] = dm_model_s1.get_avg_trace();
        display(['trace section 1: ' num2str(trace10) , '; ' num2str(trace12)]);
        display(['Iteration: ' num2str(iter_ctr)]);
        display(['average bias (kV/cm) sec 1: ' num2str(mean(TL_model_s1.v_TL)*10)]);
        display(['v_1TL(1) (V): ' num2str(TL_model_s1.v_TL(1)*10)]);
        display(['i_in_1 (A): ' num2str(TL_model_s1.i_TL(1)*params_s1.width)]);
        display(['i_out_1 (A): ' num2str(J_tot1*params_s1.width)]);
    end


    %%%% obtain the field, field intensity and the total population at position "idx" ...
    %store fields info
    record_U(iter_ctr)= dat.U(idx);
    record_V(iter_ctr)= dat.V(idx);
    record_v_TL(iter_ctr)= TL_model_s1.v_TL(idx);
    record_i_TL(iter_ctr)= TL_model_s1.i_TL(idx);
    record_r110(iter_ctr) = dm_model_s1.r110(idx); 
    record_r330(iter_ctr) = dm_model_s1.r330(idx);
    record_r220(iter_ctr) = dm_model_s1.r220(idx);
    record_J_TL(iter_ctr) = J_TL1(idx);
%     record_rRES = cell(NLVLS-4,padsize);
%     record_popsum = zeros(padsize,1);

    dat.t = dat.t+dt;
    iter_ctr = iter_ctr + 1;

end

save(savename);
end

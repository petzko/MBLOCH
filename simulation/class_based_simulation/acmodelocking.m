function acmodelocking(simfilename,interpDataFile, savename,ploton, varargin)

%parse all input files and load the scatterin rates file !
sim_params_in = false;
simRT_in = false;
if length(varargin) > 0 

    for idx =1:2:length(varargin)
        argname = varargin{idx};
        argval = varargin{idx+1};

        if strcmp(argname,'simsettings') 
            sim_params = argval;
            sim_params_in = true;
        elseif strcmp(argname,'simrt')
            simRT = argval;
            simRT_in = true;
        else
            display([ 'Unknown input argument ' argname  'with value:'])
            argval;
        end
    end
    
end
if ~sim_params_in
    sim_params = input_parser(simfilename);
end
load(interpDataFile);


bias = sim_params.bias ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZE DEVICE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%
%phase velocity inside the medium ( in mm per picosecond ... )RT time and
%frequency
tch = sim_params.tch;
lch = sim_params.lch;
c = Constants('c',{'time',tch},{'length',lch})/sim_params.nTHz;

sim_params.c = c;
Ltot = sim_params.length;

T_R = 2*Ltot/c; f_R = 1/T_R;
sim_params.INJ = INJ; sim_params.ULL = ULL; sim_params.LLL = LLL; sim_params.DEPOP = DEPOP;
sim_params.IGNORELEVEL = -1; sim_params.NLVLS = NLVLS;



hbar = Constants('hbar',{'time',tch})/Constants('q0');
%central field frequency.
E0 =  (E_fit{ULL}(bias/10)-E_fit{LLL}(bias/10))/hbar;

sim_params.E0 = E0;

N = sim_params.N_pts;
modF = sim_params.modF; % modulation frequency
%grid size in x direction
x = linspace(0,Ltot,N)';
dx = x(2) - x(1);
dt = dx/c;


sim_params.dx = dx;
sim_params.dt = dt;

sim_params.zUL = zUL_fit(bias/10);

sim_params.zNORM = sim_params.zUL;

sim_params.Ncarriers_cm = sim_params.dN*sim_params.Ld/sim_params.Lp; %carrier density

sim_params.IDX = 1:sim_params.N_pts;
sim_params.x = x(sim_params.IDX);

record_U= 1;

TL_bias = bias/10*ones(N,1);

dm_model_ = DM_MODEL_PSEUDO_3_LVL_RWA_FP(sim_params);
dm_model_.interpolate(TL_bias(sim_params.IDX),W_fit,E_fit,AC_fit,zUL_fit);


dat.N = N; dat.c = c; dat.dx = dx;
dat.dt = dt;
dat = makeMaxwellVars(dat);

rlgc.C = 1;              %unit: pF/mm
rlgc.L = 1.6e2;          %unit: pH/mm
rlgc.R = 45*sqrt(modF);  %unit: Ohm/mm --from paper W. Maineult: Microwave modulation of terahertz quantum cascade lasers: a transmission-line approach

TL_model_s1 = TL_solver2(sim_params,rlgc);


%simulation info storage arrays -> preallocate
idx = round(N/2); iter_ctr = 1;
if ~simRT_in 
    dat.simRT = 400; 
else
    dat.simRT = simRT; 
end

tEnd = dat.simRT*T_R; % end time in tps

%simulation info storage arrays -> preallocate
recordingduration = tEnd; % how many ps should we record the pulse for
iterperrecord = 1; recordingiter  = round(recordingduration/iterperrecord/dt);
padsize = double(recordingiter-length(record_U));

record_U= zeros(padsize,1);
record_V = zeros(padsize,1);
record_v_TL =zeros(padsize,1);
record_i_TL = zeros(padsize,1);
record_time = zeros(padsize,1);
record_J_TL =zeros(padsize,1);
record_r110 = zeros(padsize,1);
record_r330 = zeros(padsize,1);
record_r220 = zeros(padsize,1);


dat.t = 0;
P = zeros(dat.N,1); P_t = zeros(dat.N,1); M = P; M_t = P_t; losses = P_t;
f_display = 500;

suffix ='_';

while( dat.t< tEnd)
    
    dm_model_.propagate(dat.U,dat.V,dt);
    
    [P1,P1_t,M1,M1_t,losses1] = dm_model_.get_polarization_and_losses();
    P(sim_params.IDX) = P1; P_t(sim_params.IDX) = P1_t;
    M(sim_params.IDX) = M1; M_t(sim_params.IDX) = M1_t;
    losses(sim_params.IDX) = losses1;
    
    
    J_TL1 = dm_model_.get_current_density(sim_params);
    J_tot1 = trapz(x(1:sim_params.N_pts),J_TL1);
    
    dat = stepWave(dat,P,P_t,M,M_t,losses);
    
    %   Transmission line equations
    if iter_ctr >= 2000
        if iter_ctr == 2000        %Set initial current distribution
           J_TL0 = J_TL1;
        end
        TL_model_s1.propagate2(J_TL1,J_TL0,dat.t-dt*1999)
    end
    J_TL0 = J_TL1;
    dm_model_.update_state();
    dm_model_.interpolate(TL_model_s1.v_TL,W_fit,E_fit,AC_fit,zUL_fit);
 
    if(iter_ctr >= 2000 && mod(iter_ctr,f_display) == 0)
        clc;
        trace10 = dm_model_.get_avg_trace();
        if isnan(trace10)
            suffix = '_NAN'
            break
        end
        display(['trace section 1: ' num2str(trace10) ]);
        display(['Iteration: ' num2str(iter_ctr) '/' num2str(round(tEnd/dat.dt))]);
        display(['average bias (kV/cm) sec 1: ' num2str(mean(TL_model_s1.v_TL)*10)]);
        display(['v_1TL(1) (V): ' num2str(TL_model_s1.v_TL(1)*10)]);
        display(['i_ac (A): ' num2str(TL_model_s1.i_TL(1)*sim_params.width)]);
        display(['i_out (A): ' num2str(J_tot1*sim_params.width)]);
         
        if ploton
            subplot(2,1,1);
            plot(x,abs(dat.U).^2+abs(dat.V).^2);
            xlabel('x (mm)'); ylabel('intensity (ps^{-2}');
            subplot(2,1,2);
            plotyy(x,TL_model_s1.v_TL*10,x,[TL_model_s1.i_TL,J_TL1]);
            xlabel('x (mm)'); legend('v (kV/cm)','i (A)','J (A/mm^2)');
            getframe;
            
        end
    end
    
    
    dat.t = dat.t+dt;
    iter_ctr = iter_ctr + 1;
    
    
    %%%% obtain the field, field intensity and the total population at position "idx" ...
    %store fields info
    record_U(iter_ctr)= dat.U(idx);
    record_V(iter_ctr)= dat.V(idx);
    record_v_TL(iter_ctr)= TL_model_s1.v_TL(idx)*sim_params.height*1e3;%unit: V
    record_i_TL(iter_ctr)= TL_model_s1.i_TL(idx)*sim_params.width;     %unit: A
    record_time(iter_ctr)= dat.t;
    record_J_TL(iter_ctr) = J_TL1(idx);
    record_r110(iter_ctr) = dm_model_.r110(idx);
    record_r330(iter_ctr) = dm_model_.r330(idx);
    record_r220(iter_ctr) = dm_model_.r220(idx);
    
end



save([savename suffix]);
end

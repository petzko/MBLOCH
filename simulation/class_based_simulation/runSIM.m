clear;clc;close all;
scenarioFile = 'TL_ALL.set';
ratesfile = 'Z:\workspaces\PeTz\REPOS\MBLOCH\simulation\TL SIMS FINAL\sim_files\QCL183S\data\fitted_data1.mat';

%parse all input files and load the scatterin rates file !
params = input_parser(scenarioFile);
load(ratesfile);
    
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZE DEVICE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%
%phase velocity inside the medium ( in mm per picosecond ... )RT time and
%frequency
params.c = Constants('c',{'time',params.tch},{'length',params.lch})/params.nTHz; params.T_R = 2*params.Ltot/params.c; params.f_R = 1/params.T_R;
params.INJ = INJ; params.ULL = ULL; 
params.LLL = LLL; params.DEPOP = DEPOP;
params.NLVLS = NLVLS;
params.N_rest = params.NLVLS-5;

if exist('IGNORELEVEL')
    params.IGNORELEVEL = IGNORELEVEL;
else
    params.IGNORELEVEL = -1;
end

hbar = Constants('hbar',{'time',params.tch})/Constants('q0');
%central field frequency.
params.E0 =  (E_fit{params.ULL}(params.bias/10)-E_fit{params.LLL}(params.bias/10))/hbar;
%cavity loss l_0 in (cm^-1) --> l_0*100 in (m^-1) --> 1 mm^-1
params.linear_loss = params.linear_loss*100/(1/params.lch);

params.N = 5000;
params.N_pts = params.N;

%grid size in x direction
params.x = linspace(0,params.Ltot,params.N)';
params.T_R = 2*params.Ltot/params.c;
params.f_R = 1./params.T_R;

params.dx = params.x(2) - params.x(1); 
params.dt = params.dx/params.c;
params.zUL = zUL_fit(params.bias/10);
params.zNORM = params.zUL;

params.Ncarriers_cm = params.dN*params.Ld/params.Lp; %carrier density
params.IDX = 1:params.N;


record_U= 1; record_V = 1;

record_r110 = 1; record_r330 = 1; record_r220 = 1;
record_rRES = cell(params.NLVLS-4,1);

record_popsum = 1; record_v_TL =1;
record_i_TL = 1; record_J_TL =1;


v_TL = params.bias/10*ones(params.N,1);
dm_model = DM_MODEL_3_LVL_RWA_FP(params);
dm_model.interpolate(v_TL,W_fit,E_fit,AC_fit,zUL_fit);
dat = 1;
dat = makeMaxwellVars(params,dat);
dat = makeTransLineVarsx2(params,dat);

%%%% specify some of the mainloop control parameters %%%%
idx = 1; ctr = 1; iter_ctr = 1;

iter_per_rt = round(params.T_R/params.dt); 
params.simRT = 100; tEnd = params.simRT*params.T_R; % end time in tps
interpCtr = 500; %set how often to interpolate the energies, scattering rates, dipole elements etc.
checkptIter = 100000; 

%simulation info storage arrays -> preallocate
recordingduration = tEnd; % how many ps should we record the pulse for
iterperrecord = 1; recordingiter  = round(recordingduration/iterperrecord/params.dt);
padsize = double(recordingiter-length(record_U));

%preallocate memory for optical field, curr density voltage wave, etc.
%strorage...
record_U= padarray(record_U,padsize,'post'); record_V = padarray(record_V,padsize,'post');
record_v_TL = padarray(record_v_TL,padsize,'post');
record_i_TL = padarray(record_i_TL,padsize,'post');
record_J_TL = padarray(record_J_TL,padsize,'post');
%store population info
record_r110 = padarray(record_r110,padsize,'post');
record_r330 = padarray(record_r330,padsize,'post');
record_r220 = padarray(record_r220,padsize,'post');
for i = 1:params.N_rest
    record_rRES{i} = padarray(record_rRES{i},padsize,'post');
end
record_popsum = padarray(record_popsum,padsize,'post');

info.settings = params;
info.cavity = 'FP-QCL183S';
info.Ltot = params.Ltot;
info.N = params.N;
info.SIMTYPE = ['TL with (i0,v0) = (' num2str(params.current*10) ',' num2str(params.bias) ') in units (A/cm,kV/cm)'];

t = 0;
while( t< tEnd)
    
    
    
    dm_model.propagate(dat.U,dat.V,params.dt);
    [P,P_t,M,M_t,losses] = dm_model.get_polarization_and_losses();
    J_TL = dm_model.get_current_density(params);
    
    dat = stepWave(params,dat,P,P_t,M,M_t,losses);
    dat = stepTransLinex2(params,dat,J_TL);
    dm_model.update_state();
    
    if mod(iter_ctr,interpCtr) == 0
        dm_model.interpolate(dat.v_TL,W_fit,E_fit,AC_fit,zUL_fit);
    end
    

    %%%%% end of setting up the TL PARAMS %%%%%
    
    if(mod(iter_ctr+1,checkptIter) == 0 )
        checkptname = ['CHCKPT_' params.name '_N_TRANSMISSION_LINE_' num2str(params.N_pts) '_FP'];
        save(checkptname);
    end
    
    %%plot some of the results if neeed ariseth :D
    if(mod(iter_ctr,500) == 0)
        clc;
        info.iter_ctr = iter_ctr;
        info.RT = t/params.T_R;
        intensity = dat.U.*conj(dat.U) + dat.V.*conj(dat.V) ;
        info.maxInt  =  max(intensity);
        printINFO(info);
        display(['Vs: ' num2str(dat.Vs(t))]);
        display(['v_0: ' num2str(dat.v_TL(1))]);
        display(['i-in: ' num2str(dat.i_TL(1))]);
        display(['i-out_1: ' num2str(trapz(params.x,J_TL))]);
        display(['trace = ' num2str(popsum)])
        
        subplot(3,1,1)
        plot(params.x,[real(dat.U),real(dat.V)]);
        subplot(3,1,2)
        ax = plotyy(params.x(1:params.N_pts-1),[dat.v_TL(1:params.N_pts-1)*10], ...
        params.x(1:params.N_pts-1),dat.i_TL(1:params.N_pts-1));
        title(info.SIMTYPE);
        subplot(3,1,3)
        %plots the populations and the current density
        plotyy(params.x,[dm_model.r110,dm_model.r330,dm_model.r220],params.x,J_TL);
        legend('\rho_{11}','\rho_{33}','\rho_{22}','J');
        
        getframe;
    end
    
    
    %%%% obtain the field, field intensity and the total population at position "idx" ...
    %store fields info
    record_U(ctr)= dat.U(idx);  record_V(ctr)= dat.V(idx);
    record_v_TL(ctr) = dat.v_TL(idx);
    record_i_TL(ctr) = dat.i_TL(idx);

    record_J_TL(ctr) = J_TL(idx);

    record_r110(ctr) = dm_model.r110(idx);
    record_r220(ctr) = dm_model.r220(idx);
    record_r330(ctr) = dm_model.r330(idx);
    popsum = dm_model.r110(idx)+dm_model.r220(idx)+dm_model.r330(idx) ;
    for i = 1:dm_model.N_rest
        record_rRES{i}(ctr) = dm_model.rRES{i}(idx);
        popsum = popsum + dm_model.rRES{i}(idx);
    end
    record_popsum(ctr) =  popsum;
    ctr = ctr+1;
    
    t = t+params.dt;
    iter_ctr = iter_ctr + 1;
    
end
savename = [params.name '_' params.scenario '_N_' num2str(params.N) '_FP_' num2str(params.simRT) ];
save(savename);


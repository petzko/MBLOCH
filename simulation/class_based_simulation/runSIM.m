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
params_s1.NLVLS = NLVLS; params_s1.N_rest = params_s1.NLVLS-4;

params_s2.INJ = INJ; params_s2.ULL = ULL; params_s2.LLL = LLL; params_s2.DEPOP = DEPOP;
params_s2.NLVLS = NLVLS; params_s2.N_rest = params_s2.NLVLS-4;


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

dm_model_s1 = DM_MODEL_3_LVL_RWA_FP(params_s1);
dm_model_s1.interpolate(TL_bias(params_s1.IDX),W_fit,E_fit,AC_fit,zUL_fit);

dm_model_s2 = DM_MODEL_3_LVL_RWA_FP(params_s2);
dm_model_s2.interpolate(TL_bias(params_s2.IDX),W_fit,E_fit,AC_fit,zUL_fit);

dat.N = N; dat.c = c; dat.dx = dx;
dat.dt = dt;
dat = makeMaxwellVars(dat);

rlgc1.C = 1.3;    %unit: pF/mm
rlgc1.L = 1.6e2;  %unit: pH/mm
rlgc2.C = 1.3;    %unit: pF/mm
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

t = 0;
P = zeros(dat.N,1); P_t = zeros(dat.N,1); M = P; M_t = P_t; losses = P_t;
interpCtr = 1; %set how often to interpolate the energies, scattering rates, dipole elements etc.
checkptIter = 1040000;% 1039800; %100000
f_plot = 1000;
f_display = 100;


while( t< tEnd)
    
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

    dat = stepWave(dat,P,P_t,M,M_t,losses);

%   Transmission line equations
    if iter_ctr > 1000
        TL_model_s1.propagate_1(J_TL1);
        TL_model_s2.propagate_2(J_tot1,J_TL2);
    end


    dm_model_s1.update_state();
    dm_model_s2.update_state();
    
     if mod(iter_ctr,interpCtr) == 0 %interpCtr
            MM1 = max(TL_model_s1.v_TL); NN1 = min(TL_model_s1.v_TL);%9kV/cm--12kV/cm
            if MM1 < 1.2 && NN1 > 0.9
            VV_TL1 = TL_model_s1.v_TL;
            end
             
            MM2 = max(TL_model_s2.v_TL); NN2 = min(TL_model_s2.v_TL);%9kV/cm--12kV/cm
            if MM2 < 1.2 && NN2 > 0.9
            VV_TL2 = TL_model_s2.v_TL;
            end

        dm_model_s1.interpolate(VV_TL1,W_fit,E_fit,AC_fit,zUL_fit);
        dm_model_s2.interpolate(VV_TL2,W_fit,E_fit,AC_fit,zUL_fit);
    end

    if(mod(iter_ctr,f_display) == 0)
        clc;
        trace1 = dm_model_s1.get_avg_trace();
        trace2 = dm_model_s2.get_avg_trace();
        display(['trace section 1: ' num2str(trace1)]);
        display(['trace section 2: ' num2str(trace2)]);
        display(['Iteration: ' num2str(iter_ctr)]);
        display(['average bias (kV/cm) sec 1: ' num2str(mean(TL_model_s1.v_TL)*10)]); 
        display(['average bias (kV/cm) sec 2: ' num2str(mean(TL_model_s2.v_TL)*10)]); 
        display(['v_1TL(1) (V): ' num2str(TL_model_s1.v_TL(1)*10)]);
% % % %         display(['i_1TL(1) (A): ' num2str(i1new)]);
        display(['v_2TL(1) (V): ' num2str(TL_model_s2.v_TL(1)*10)]);
% % % %         display(['i_2TL(1) (A): ' num2str(i2new)]);

% % % %         i_out = trapz(x(1:params_s1.N_pts),J_TL1)*params_s1.width;
        display(['i_out (A): ' num2str(J_tot1*params_s1.width)]);
%     end
% 
%     
%     %%plot some of the results if neeed ariseth :D
%     if(mod(iter_ctr,f_plot) == 0)

        subplot(4,1,1);
        plot(x,abs(dat.U).^2,x,abs(dat.V).^2);
        title(['t=',num2str(t),' ps','        |E_z|^2','        Iteration= ',num2str(iter_ctr)]);
        
        subplot(4,1,2); 
%         J_TL(dm_model_s1.IDX) = J_TL1; J_TL(dm_model_s2.IDX) = J_TL2;
        J_TL(1:params_s1.N_pts) = J_TL1; J_TL(params_s1.N_pts+1:params_s1.N_pts+params_s2.N_pts) = J_TL2;        
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
%         frame = getframe(1);
%         GIF{iter_ctr/f_plot} = frame2im(frame);
    end
    
    
    
    %%%%% end of setting up the TL PARAMS %%%%%
%     if(mod(iter_ctr+1,checkptIter) == 0 )
%         checkptname = ['CHCKPT_' params.name '_N_TRANSMISSION_LINE_' num2str(params.N_pts) '_FP'];
%         save(checkptname);
%     end
    
    
    %%%% obtain the field, field intensity and the total population at position "idx" ...
    %store fields info
    record_U(iter_ctr)= dat.U(idx);  record_V(ctr)= dat.V(idx);
    
    t = t+dt;
    iter_ctr = iter_ctr + 1;
    
end



% % Creat gif image
% filename = 'QCL dynamics.gif'; % Specify the output file name
%         for idx = 1:iter_ctr/f_plot
%             [A,map] = rgb2ind(GIF{idx},256);
%             if idx == 1
%             imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.5);
%             else
%             imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.5);
%             end
%         end

% Fourier transformation 
mydft(record_U,dt);

simname = 'sim-test';
        
savename = [simname];
save(savename);


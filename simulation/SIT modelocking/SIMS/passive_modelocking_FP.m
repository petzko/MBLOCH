function passive_modelocking(params_gain,params_abs,pump_strength,savename,model)

%%
assert(params_gain.tch == params_abs.tch,'Characteristic time mismatch');
assert(params_gain.lch == params_abs.lch,'Characteristic length mismatch');
% resonance frequency in the absorber and the gain;
E0_1 = (params_abs.E_UL-params_abs.E_LL);
E0_2 = (params_gain.E_UL-params_gain.E_LL);

% the assumed central frequency of the electric field
params_gain.E0 = 1/2*(E0_1+E0_2); % (angular freq)
params_abs.E0 = 1/2*(E0_1+E0_2); %  (angular freq)
% total number of grid points
N = 3000;
% total length (sum/max?)
Ltot = params_abs.L+params_gain.L;

%%
%{
now here we calculate the position of the gain and the abosrber along the
cavity; We can do things such as:
[gain;abosrber]
[absorber;gain]
or even a grating like structure
[gain;absorber;gain;absorber;gain ... ].
%}

% fraction of grid points in the gain section and the abs section
%{
NOW HERE IS where we decide the ordering [gain;absorber] or [absorger;gain]
or even something more exotic; Simply set the IDX vector of the
corresponding material telling us for which grid points is the density
matrix valid; E.g. for [abs;gain] we set:
 params_abs.IDX = 1:params_abs.N_pts;
    and
 params_gain.IDX = params_abs.N_pts+1:N;
%}


N_gain = round(params_gain.L/Ltot*N);  params_gain.N_pts = N_gain;
N_abs = N-N_gain; params_abs.N_pts = N_abs;
IDX = [ceil(params_abs.N_pts/2), params_gain.N_pts,floor(params_abs.N_pts/2)];
IDX = cumsum(IDX); 

params_abs.IDX = [1:IDX(1),IDX(2)+1:IDX(3)];
params_gain.IDX = [IDX(1)+1:IDX(2)];


%total grid vector!
x = linspace(0,Ltot,N);
dx = x(2)-x(1);


% we need the background refractive index to set the time-step
nThz = 1/2*(params_gain.nTHz+params_abs.nTHz);
% phase velocity at the envelope central frequency;
c = Constants('c',{'time',params_abs.tch},{'length',params_abs.lch})/nThz;
% The condition dt = dx/c is imposed by our numerical scheme -> to achieve second order
% accuracy we need dt = dx/c;
dt = dx/c;

% round trip time
T_RT = 2*Ltot/c;
% field envelope vector!
F = 0.001*rand(N,1);
B = 0.001*rand(N,1);
F_solver = RNFDSolver(N,dx,1,c,F);
B_solver = RNFDSolver(N,dx,-1,c,B);

% with respect to which dipole transition do you normalize the field !?!?!
zNORM = params_gain.zUL;
% propagate this info
params_abs.zNORM = zNORM;
params_gain.zNORM = zNORM;


%%
%{
Now here we initialize the density matrix solvers. Those are in turn coupled to
the optical field via the electric dipole interaction (zUL) and the
polarization term. The abbreviation tells us that DM_MODEL_2_LVL_RWA_RING
tells us that this class models a 2 level denstiy matrix system, within the
rotating wave approximation on a ring cavity laser. We assume one two level
system per grid point. Here we can also use different models for the laser.
For example use a three level model, which coherently includes the resonant
tunneling into the dynamics;
%}

params_gain = DM_MODEL_2_LVL_RWA_FP.calc_threshold_inv(params_gain,params_abs,pump_strength);
d_th = params_gain.d_th;
gain_model = DM_MODEL_2_LVL_RWA_FP(params_gain);
WG_sat = 1/(gain_model.T1*gain_model.T2);


%or this for a two level system model for the gain
%absorber model
abs_model = DM_MODEL_2_LVL_RWA_FP(params_abs);

% these vectors are initialized for plotting purposes;
ruu = zeros(N,1); rll = zeros(N,1); eta_ = rll;
%%
%{
This is the main loop
%}
t = 1; N_t = 2000000;

%{
Initialize variable storage arrays
%}

%preallocate memory for optical field, curr density voltage wave, etc.
%strorage...
% field variable at the end of the gain section and the absorber section;
record_F_g= zeros(N_t,1); record_F_a = zeros(N_t,1);
record_B_g= zeros(N_t,1); record_B_a = zeros(N_t,1);

%store population info

% gain and absorver pop densities
record_r22g = zeros(N_t,1); record_r11g = zeros(N_t,1);
record_r22a = zeros(N_t,1); record_r11a = zeros(N_t,1);
record_n21g = zeros(N_t,1); record_n21a = zeros(N_t,1);

record_INV_G = zeros(N_t,1);
record_INV_A = zeros(N_t,1);
record_TGR = zeros(N_t,1);


info.simname = 'FAST SA PML SIMULATIONS.';
info.cavity = 'RING';
info.simRT = N_t*dt/(T_RT);
info.Ltot = Ltot;
info.N = N;
info.p = pump_strength;
info.d_th = d_th;
info.SIMTYPE = ['Abs length: ' num2str(params_abs.L) 'mm; Gain length: ' num2str(params_gain.L) ' mm'];
d_th_vector = d_th * ones(length(x),1);

iter_per_rt = round(T_RT/dt);

last_5_rts_ids = N_t-5.2*iter_per_rt;
skip = 5;
F_data_vs_x = zeros(length(x(1:skip:end)),2*iter_per_rt);
B_data_vs_x = zeros(length(x(1:skip:end)),2*iter_per_rt);

inversion_data_gain_vs_x = zeros(length(gain_model.IDX(1:skip:end)),2*iter_per_rt);
inversion_data_abs_vs_x = zeros(length(abs_model.IDX(1:skip:end)),2*iter_per_rt);

additional_ctr = 1;

while t <= N_t
    
    ruu(gain_model.IDX) = gain_model.rho_u_DC;
    rll(gain_model.IDX) = gain_model.rho_l_DC;
    
    ruu(abs_model.IDX) = abs_model.rho_u_DC;
    rll(abs_model.IDX) = abs_model.rho_l_DC;
    % store info!
    delta_0 = gain_model.rho_u_DC - gain_model.rho_l_DC;
    if strcmp(model,'2LVL')
        T_GR = gain_model.T1*log((pump_strength*params_gain.d_th-delta_0)...
            ./(params_gain.d_th*(pump_strength-1)));
        record_TGR(t) = T_GR(end);
    end
    
    %store fields info
    record_F_g(t) = F(gain_model.IDX(end));
    record_B_g(t) = B(gain_model.IDX(end));
    
    record_r22g(t) = gain_model.rho_u_DC(end);
    record_r11g(t) = gain_model.rho_l_DC(end);
    
    record_INV_G(t) = dx*trapz(gain_model.rho_u_DC-gain_model.rho_l_DC);
    record_INV_A(t) = dx*trapz(abs_model.rho_u_DC-abs_model.rho_l_DC);
    record_F_a(t)= F(abs_model.IDX(end));
    record_B_a(t)= B(abs_model.IDX(end));
    
    record_r22a(t) = abs_model.rho_u_DC(end);
    record_r11a(t) = abs_model.rho_l_DC(end);
    
    % record at every grid point for the last two round trips
    if t >= last_5_rts_ids
        F_data_vs_x(:,additional_ctr) = F(1:skip:end);
        B_data_vs_x(:,additional_ctr) = B(1:skip:end);
        inversion_data_gain_vs_x(:,additional_ctr) = gain_model.rho_u_DC(1:skip:end)-gain_model.rho_l_DC(1:skip:end);
        inversion_data_abs_vs_x(:,additional_ctr) = abs_model.rho_u_DC(1:skip:end)-abs_model.rho_l_DC(1:skip:end);
        
        additional_ctr = additional_ctr+1;
    end
    
    
    %     % timestep the von neumann equation for the gain/absorber density matrices;
    %     gain_model.propagate(F(gain_model.IDX),B(gain_model.IDX),dt); %  select the only field values that overlap with the gain
    %     abs_model.propagate(F(abs_model.IDX),B(abs_model.IDX),dt);
    %
    %     % polarization, it's first order time derivative and the loss vector;
    %     % note that we need to reinit every time in case we want to test a
    %     % situation when the gian and the absorber are overlapping;
    P = zeros(N,1); P_t = zeros(N,1);
    losses = zeros(N,1);
    M = zeros(N,1); M_t = zeros(N,1);
    
    gain_model.propagate(F(gain_model.IDX),B(gain_model.IDX),dt);
    abs_model.propagate(F(abs_model.IDX),B(abs_model.IDX),dt);
    
    [Gdummy1,Gdummy2,Gdummy3,Gdummy4,Gdummy5] = gain_model.get_polarization_and_losses();
    P(gain_model.IDX) = Gdummy1;
    P_t(gain_model.IDX) = Gdummy2;
    
    M(gain_model.IDX) = Gdummy3;
    M_t(gain_model.IDX) =Gdummy4;
    
    losses(gain_model.IDX) =  Gdummy5; %typing mistake dummy5 instead of dummy3
    
    [dummy1,dummy2,dummy3,dummy4,dummy5] = abs_model.get_polarization_and_losses();
    P(abs_model.IDX) = P(abs_model.IDX) + dummy1;
    P_t(abs_model.IDX) = P(abs_model.IDX) + dummy2;
    
    M(abs_model.IDX) = M(abs_model.IDX) + dummy3;
    M_t(abs_model.IDX) = M_t(abs_model.IDX)  + dummy4;
    
    losses(abs_model.IDX) = dummy5;
    
    F = F_solver.make_step(-1i*c*P,-1i*c*P_t,-c*losses,dt);
    B = B_solver.make_step(-1i*c*M,-1i*c*M_t,-c*losses,dt);
    
    %%% The solver make_step(obj,F,F_t,k,dt) assumes an equation of the form
    %%% D_t(E) = -/+ D_x(E) + F(x,t) + k.*E
    F = F_solver.set_bdry(B(1),'no');
    B = B_solver.set_bdry('no',F(end));
    gain_model.update_state(); 
    abs_model.update_state(); 
    
    % now do some plotting
    if mod(t,1000) == 0
        clc;
        info.iter_ctr = t;
        info.RT = t*dt/T_RT;
        intensity_F = F.*conj(F);
        intensity_B = B.*conj(B);
        
        info.maxInt  =  max(intensity_F);
        printINFO(info);
        
        subplot(311);
        plot(x,intensity_F,'-b',x,intensity_B);
        xlabel('x (mm)'); ylabel('|F|^2,|B|^2 (arb. u.)');
        
        
        xlabel('x (mm)'); ylabel('pop. inversion(s)');
        
        if strcmp(model,'2LVL')
            subplot(312);
           
            plot(x,ruu-rll);
            subplot(313);
            delta_0_predicted= pump_strength*params_gain.d_th./(1+abs(F(gain_model.IDX)).^2/WG_sat);
            ax =  plotyy(x(gain_model.IDX),[delta_0,delta_0_predicted], ...
                x(gain_model.IDX),[real(T_GR),T_RT*ones(gain_model.N_pts,1)]);
            legend('\Delta_0(t)','\Delta_{min}','\tau_{gr}','T_{rt}');
            ax(1).YLim = [0,1];
            ax(2).YLim = [0,T_RT+4];
        else
            subplot(312);
            
            plot(x(gain_model.IDX),gain_model.rho_u-gain_model.rho_l,x(abs_model.IDX),abs_model.rho_u-abs_model.rho_l,x,d_th_vector);
            legend('gain inversion','abs inversion','\Delta_{th}');
        end
        getframe;
        
    end
    
    
    t = t+1;
    
end
save(savename);

end




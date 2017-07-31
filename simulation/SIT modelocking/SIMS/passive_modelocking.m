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

overlap = false
if ~overlap
    
    N_gain = round(params_gain.L/Ltot*N);  params_gain.N_pts = N_gain;
    N_abs = N-N_gain; params_abs.N_pts = N_abs;
    
    params_gain.IDX = 1:params_gain.N_pts;
    params_abs.IDX = params_gain.N_pts+1:N;
else
    params_gain.L = Ltot;
    params_abs.L = Ltot;
    
    params_gain.N_pts = N;
    params_abs.N_pts = N;
    params_gain.IDX = 1:N;
    params_abs.IDX = 1:N;
    
end

%total grid vector!
x = linspace(0,Ltot,N);
dx = x(2)-x(1);
% field envelope vector!
F = zeros(N,1);

% we need the background refractive index to set the time-step
nThz = 1/2*(params_gain.nTHz+params_abs.nTHz);
% phase velocity at the envelope central frequency;
c = Constants('c',{'time',params_abs.tch},{'length',params_abs.lch})/nThz;
% The condition dt = dx/c is imposed by our numerical scheme -> to achieve second order
% accuracy we need dt = dx/c;
dt = dx/c;

% round trip time
T_RT = Ltot/c;

F_solver = RNFDSolver(N,dx,1,c,F);
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

% use this for a 3lvl system model for the gain with RT included
if strcmp(model,'RT')
    
    gain_model = DM_MODEL_3_LVL_RWA_RING_v2(params_gain);
    d_th = gain_model.calc_threshold_inv(params_gain,params_abs);
    
    inversion = gain_model.steady_state(gain_model.ULL)-gain_model.steady_state(gain_model.LLL);
    pump_strength = inversion/d_th;
    
    %or this for a two level system model for the gain
elseif strcmp(model,'MLVL')
    gain_model = DM_MODEL_MLVL_RWA_RING(params_gain);
    d_th = gain_model.calc_threshold_inv(params_gain,params_abs);
    
    inversion = gain_model.steady_state(gain_model.ULL)-gain_model.steady_state(gain_model.LLL);
    pump_strength = inversion/d_th;
else
    params_gain = DM_MODEL_2_LVL_RWA_RING.calc_threshold_inv(params_gain,params_abs,pump_strength);
    d_th = params_gain.d_th;
    gain_model = DM_MODEL_2_LVL_RWA_RING(params_gain);
    WG_sat = 1/(gain_model.T1*gain_model.T2);
   
end

%or this for a two level system model for the gain
%absorber model
abs_model = DM_MODEL_2_LVL_RWA_RING(params_abs);

% these vectors are initialized for plotting purposes;
ruu = zeros(N,1); rll = zeros(N,1); eta_ = rll;
%%
%{
This is the main loop
%}
t = 1; N_t = 1000000;

%{
Initialize variable storage arrays
%}

%preallocate memory for optical field, curr density voltage wave, etc.
%strorage...
% field variable at the end of the gain section and the absorber section;
record_U_g= zeros(N_t,1); record_U_a = zeros(N_t,1);
%store population info

% gain and absorver pop densities
record_r22g = zeros(N_t,1); record_r11g = zeros(N_t,1);
record_r22a = zeros(N_t,1); record_r11a = zeros(N_t,1);
record_n21g = zeros(N_t,1); record_n21a = zeros(N_t,1);

record_INV_G = zeros(N_t,1); record_INV_A = zeros(N_t,1);
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
pulse_data_vs_x = zeros(length(x(1:skip:end)),2*iter_per_rt);
inversion_data_gain_vs_x = zeros(length(gain_model.IDX(1:skip:end)),2*iter_per_rt);
inversion_data_abs_vs_x = zeros(length(abs_model.IDX(1:skip:end)),2*iter_per_rt);

additional_ctr = 1;

while t <= N_t
    
    
    % store info!
    delta_0 = gain_model.rho_u - gain_model.rho_l;
    if strcmp(model,'2LVL')
        T_GR = gain_model.T1*log((pump_strength*params_gain.d_th-delta_0)...
            ./(params_gain.d_th*(pump_strength-1)));
        record_TGR(t) = T_GR(end);
    end
    
    %store fields info
    record_U_g(t) = F(gain_model.IDX(end));
    record_r22g(t) = gain_model.rho_u(end);
    record_r11g(t) = gain_model.rho_l(end);
    record_n21g(t) = gain_model.eta_ul(end);
    
    record_INV_G(t) = dx*trapz(gain_model.rho_u-gain_model.rho_l);
    record_INV_A(t) = dx*trapz(abs_model.rho_u-abs_model.rho_l);
    record_U_a(t)= F(abs_model.IDX(end));
    record_r22a(t) = abs_model.rho_u(end);
    record_r11a(t) = abs_model.rho_l(end);
    record_n21a(t) = abs_model.eta_ul(end);
    
    % record at every grid point for the last two round trips
    if t >= last_5_rts_ids
        pulse_data_vs_x(:,additional_ctr) = F(1:skip:end);
        
        inversion_data_gain_vs_x(:,additional_ctr) = gain_model.rho_u(1:skip:end)-gain_model.rho_l(1:skip:end);
        inversion_data_abs_vs_x(:,additional_ctr) = abs_model.rho_u(1:skip:end)-abs_model.rho_l(1:skip:end);
        
        additional_ctr = additional_ctr+1;
    end
    
    
    % timestep the von neumann equation for the gain/absorber density matrices;
    gain_model.propagate(F(gain_model.IDX),dt); %  select the only field values that overlap with the gain
    abs_model.propagate(F(abs_model.IDX),dt);
    
    % polarization, it's first order time derivative and the loss vector;
    % note that we need to reinit every time in case we want to test a
    % situation when the gian and the absorber are overlapping;
    P = zeros(N,1); P_t = zeros(N,1); losses = zeros(N,1);
    
    [dummy1,dummy2,dummy3] = abs_model.get_polarization_and_losses();

    P(abs_model.IDX) = dummy1;
    P_t(abs_model.IDX) = dummy2;
    losses(abs_model.IDX) = dummy3;
    
    [dummy1,dummy2,dummy3] = gain_model.get_polarization_and_losses();

    P(gain_model.IDX) = P(gain_model.IDX) + dummy1; % we need to add those vectors in case P is non zero.
    P_t(gain_model.IDX) = P_t(gain_model.IDX) +dummy2;
    losses(gain_model.IDX) = losses(gain_model.IDX) + dummy3;
%     
    F = F_solver.make_step(-1i*c*P,-1i*c*P_t,-c*losses,dt);
    F = F_solver.set_bdry(F(end),'no');
    
    gain_model.update_state();
    abs_model.update_state();
    
    
    % now do some plotting
    if mod(t,1000) == 0
        clc;
        info.iter_ctr = t;
        info.RT = t*dt/T_RT;
        intensity = F.*conj(F) ;
        info.maxInt  =  max(intensity);
        printINFO(info);
        
        subplot(311);
        plot(x,intensity,'-b','Linewidth',2.0);
        xlabel('x (mm)'); ylabel('|F|^2 (arb. u.)');
        
        
        xlabel('x (mm)'); ylabel('pop. inversion(s)');
        
        if strcmp(model,'2LVL')
            subplot(312);
            positional_gain = gain_model.calculate_positional_gain(3.6)*1e-2 ;
            positional_loss = abs_model.calculate_positional_gain(3.6)*1e-2*-1;
            plot(x(gain_model.IDX),positional_gain,x(abs_model.IDX),positional_loss);
            legend('gain ','loss');
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




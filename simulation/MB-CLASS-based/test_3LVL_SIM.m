clear; clc;close all;

params_gain = input_parser('RT_GAIN_INPUTS.sim');
params_sl = input_parser('RT_SLOWLIGHT.sim');

assert(params_gain.tch == params_sl.tch,'Characteristic time mismatch');
assert(params_gain.lch == params_sl.lch,'Characteristic length mismatch');
tch = params_sl.tch; 
lch = params_sl.lch; 
Ex_ = params_sl.Ex_; 
Grnd_ = params_sl.Grnd_;
Spin_ = params_sl.Spin_;

NLVLS = 3;
params_sl.H = reshape(params_sl.H,NLVLS,NLVLS).';
params_sl.H(Spin_,Spin_) = params_sl.H(Ex_,Ex_);

params_sl.W = reshape(params_sl.W,NLVLS,NLVLS).';

E0_sd = (params_sl.H(Ex_,Ex_)-params_sl.H(Grnd_,Grnd_));
E0_gain = (params_gain.E_UL-params_gain.E_LL);

params_gain.cavity = 'RING';
params_sl.cavity = 'RING';

params_gain.E0 = 1/2*(E0_sd+E0_gain); 
params_sl.E0 = 1/2*(E0_sd+E0_gain); 
N = 2000; 
glob_idx = 1:N;

% number of regions
NR = 100; Ltot = params_sl.L;
x = linspace(0,Ltot,N)';
dx = x(2)-x(1); 

% total length of slowing regions
params_sl.L =  Ltot*0.5;
dW = params_sl.L/NR; 

% length of free space
L_GAIN = (Ltot-params_sl.L);
params_gain.L = L_GAIN;

%width of a single Free space region
dF = params_gain.L/(NR+1);
% free-space slow-down free-space slow-down free-space
core_i = zeros(length(x),NR);
core = 0*x;
for i = 1:NR
    core_i(:,i) = (x>= (i*dF+(i-1)*dW)  & x <= (i*dF+i*dW));
    core = core + core_i(:,i);
end
% interaction area !
i_core = (x>(dF) & x<=(NR*dF+NR*dW));
% slow light and gain cores 
params_sl.IDX = glob_idx(logical(core)); 
params_gain.IDX = glob_idx(logical(1-core));

N_gain = length(params_gain.IDX); 
params_gain.N_pts = N_gain;
N_sl = length(params_sl.IDX); 
params_sl.N_pts = N_sl;

F = zeros(N,1);
c = Constants('c',{'time',tch},{'length',lch})/params_sl.nTHz;
dt = dx/c; 
%%%%%%
FWHM = 30*2/sqrt(2); %psec
t0 = +50;
sig = FWHM/(2*sqrt(2*log(2))); % std. dev
A0 = 0.001;
gauss = @(t,mu,s) A0*exp(-(t-mu).^2/(2*s^2));

NBits =  8;
NUM2SEND = 125;
code = de2bi(NUM2SEND,NBits);
message = @(t) 0;

for i = 1:NBits
    message = @(t) message(t) + code(i)*gauss(t,t0+i*4*FWHM,sig)
end
%%%%%%%%%%%%
U_in = []; 
U_out = []; 


F_solver = RNFDSolver(N,dx,1,c,F);
zNORM = params_sl.zUL; % with respect to which dipole transition do you
% normalize the field !?!?! 


params_sl.zNORM = zNORM;
params_gain.zNORM = zNORM; 

gain_model = DM_MODEL_2_LVL_RWA_RING(params_gain); 
sl_model = DM_MODEL_3_LVL_RWA_RING(params_sl);
P = zeros(N,1); P_t = zeros(N,1); losses = zeros(N,1); 

t = 1; N_t = 1000000;
U_in = zeros(N_t,1);  U_out = zeros(N_t,1);

while t < N_t
    
    if mod(t,1000) == 0
        plot(x,abs(F).^2,'b','Linewidth',2.0);
%         plot(sl_model.rho_ex-sl_model.rho_grnd);
        
        xlabel('x (mm)');ylabel('Intensity (a.u.)');
        hold on;
        fobj = fill(x,(1-core)*A0^2,'g'); 
        set(fobj,'EdgeColor','none','FaceAlpha',0.3);
        hold off;
        getframe;
        clc;
        display(['gain:' num2str(gain_model.calculate_optical_gain(3.6)*1e-2)] )
    
    
    end
    U_in(t) = F(1); 
    U_out(t) = F(end);
    
    
    F_eff = F.*i_core;
    gain_model.propagate(F_eff(gain_model.IDX),dt);
    sl_model.propagate(F_eff(sl_model.IDX),dt); 
    ftr = ones(length(F),1); 
    ftr(sl_model.IDX) = -1i*c*sl_model.NORM_FACTOR_FIELD;
    ftr(gain_model.IDX) = -1i*c*gain_model.NORM_FACTOR_FIELD;
    
    [dummy1,dummy2,dummy3] = sl_model.get_polarization_and_losses();
    P(sl_model.IDX) = dummy1; 
    P_t(sl_model.IDX) = dummy2; 
    losses(sl_model.IDX) = dummy3; 
    
    [dummy1,dummy2,dummy3] = gain_model.get_polarization_and_losses();
    P(gain_model.IDX) = dummy1; 
    P_t(gain_model.IDX) = dummy2;
    losses(gain_model.IDX) = dummy3; 
    
    F = F_solver.make_step(-1i*c*P.*i_core,-1i*c*P_t.*i_core,-c*losses.*i_core,dt); 
    F = F_solver.set_bdry(message(t*dt),'no');
    
    gain_model.update_state();
    
    sl_model.update_state();
    
    F = F;
    t = t+1;
end





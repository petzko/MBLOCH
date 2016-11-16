params_gain = input_parser('GAIN_INPUTS.sim');
params_abs = input_parser('ABSORBER_INPUTS.sim');

assert(params_gain.tch == params_abs.tch,'Characteristic time mismatch');
assert(params_gain.lch == params_abs.lch,'Characteristic length mismatch');
tch = params_abs.tch; 
lch = params_abs.lch; 
E0_1 = (params_abs.E_UL-params_abs.E_LL);
E0_2 = (params_gain.E_UL-params_gain.E_LL);


params_gain.cavity = 'RING';
params_abs.cavity = 'RING';


params_gain.E0 = 1/2*(E0_1+E0_2); 
params_abs.E0 = 1/2*(E0_1+E0_2); 
N = 5000; 
Ltot = params_abs.L+params_gain.L;

N_gain = floor(params_gain.L/Ltot*N); 
params_gain.N_pts = N_gain;

N_abs = floor(params_abs.L/Ltot*N); 
params_abs.N_pts = N_abs;

N = N_abs+N_gain; 
x = linspace(0,Ltot,N); 
dx = x(2)-x(1); 
F = zeros(N,1);
c = Constants('c',{'time',tch},{'length',lch})/params_abs.nTHz;
dt = dx/c; 
tp = .1;
A0  = 1/tp;
t0 = tp*10;
F_solver = RNFDSolver(N,dx,1,c,F);
zNORM = params_gain.zUL; % with respect to which dipole transition do you
% normalize the field !?!?! 


params_abs.zNORM = zNORM;
params_gain.zNORM = zNORM; 


params_abs.IDX = 1:N_abs;
params_gain.IDX = N_abs+1:N;


gain_model = DM_MODEL_2_LVL_RWA_RING(params_gain); 
abs_model = DM_MODEL_2_LVL_RWA_RING(params_abs); 
ruu = zeros(N,1); 
rll = zeros(N,1); 
P = zeros(N,1); P_t = zeros(N,1); losses = zeros(N,1); 

t = 1; N_t = 100000;
while t < N_t
    
    if mod(t,100) == 0
        ruu(abs_model.IDX) = abs_model.rho_u;
        ruu(gain_model.IDX) = gain_model.rho_u;
        
        rll(abs_model.IDX) = abs_model.rho_l;
        rll(gain_model.IDX) = gain_model.rho_l;
        
        plotyy(x,abs(F).^2,x,ruu-rll);
        getframe;
    end
    
    gain_model.propagate(F(gain_model.IDX),dt);
    abs_model.propagate(F(abs_model.IDX),dt); 
       
    [dummy1,dummy2,dummy3] = abs_model.get_polarization_and_losses();
    P(abs_model.IDX) = dummy1; 
    P_t(abs_model.IDX) = dummy2; 
    losses(abs_model.IDX) = dummy3; 
    
    [dummy1,dummy2,dummy3] = gain_model.get_polarization_and_losses();
    P(gain_model.IDX) = dummy1; 
    P_t(gain_model.IDX) = dummy2;
    losses(gain_model.IDX) = dummy3; 
    
    F = F_solver.make_step(-1i*c*P,-1i*c*P_t,-c*losses,dt); 
    F = F_solver.set_bdry(A0*sech((t*dt-t0)/tp)+ ...
        F(end)*params_abs.R,'no');

    gain_model.update_state();
    abs_model.update_state();
    
 
    F = F;
    t = t+1;
end





powerloss_mm = powerloss*1e-3; 
sigma_g_mm = sigma_g/(1e6);
sigma_a_mm = sigma_a/(1e6);
Ng_mm = Ncarriers_g/(1e9);
Na_mm = Ncarriers_a/(1e9);


syms u v
Dl_uv = d_eq_GAIN*( (exp(-(1+2*u)*v/T1_GAIN)-1)/(1+2*u)-(exp((TR_-v)/T1_GAIN)-1))/ ... 
    (exp(-(1+2*u)*v/T1_GAIN)-exp((TR_-v)/T1_GAIN));
AVG_delta_uv = (Dl_uv-d_eq_GAIN/(1+2*u))*T1_GAIN/(1+2*u)*(1-exp(-(1+2*u)*v/T1_GAIN));
Fg_uv = sigma_g_mm*Ng_mm*2*u/TR_*AVG_delta_uv-2*powerloss_mm*u*v/TR_; 


syms w z
Dl_wz = d_eq_ABS*( (exp(-(1+2*w)*z/T1_ABS)-1)/(1+2*w)-(exp((TR_-z)/T1_ABS)-1))/ ... 
    (exp(-(1+2*w)*z/T1_ABS)-exp((TR_-z)/T1_ABS));

AVG_delta_wz = (Dl_wz-d_eq_ABS/(1+2*w))*T1_ABS/(1+2*w)*(1-exp(-(1+2*w)*z/T1_ABS));
Fa_wz = sigma_a_mm*Na_mm*2*w/TR_*AVG_delta_wz;


dFg_du = matlabFunction(diff(Fg_uv,u)); 
dFg_dv = matlabFunction(diff(Fg_uv,v)); 

dFa_dw = matlabFunction(diff(Fa_wz,w)); 
dFa_dz = matlabFunction(diff(Fa_wz,z)); 


options = optimoptions('fsolve','Display','none','PlotFcn',@optimplotfirstorderopt);
options.MaxFunEvals = 1000;
dFg = @(X) [dFg_du(X(1),X(2)) dFg_dv(X(1),X(2))];
sol_g = fsolve(dFg,[2,3],options); 

dFa = @(Y) [dFa_dw(Y(1),Y(2)),dFa_dz(Y(1),Y(2))];
sol_a = fsolve(dFa,[2,3],options); 






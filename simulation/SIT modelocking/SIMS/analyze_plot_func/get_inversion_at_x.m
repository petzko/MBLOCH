function get_inversion(x_idx,pulse_data_vs_x,T2,d_eq_GAIN)

envelope_at_X = pulse_data_vs_x(x_idx,:);
E0_sqrt_at_x = abs(envelope_at_X).^2;
E0_sqrt_max = max(E0_sqrt_at_x);
R = T2*E0_sqrt_max; % E IS NORMALIZED TO mu/hbar E not to mu/2*hbar E!!!! 


d_eq_GAIN = p*d_th;
T1_GAIN = gain_model.T1; 


inversion_at_x_idx_GAIN = inversion_data_gain_vs_x(x_idx,:);
times_at_GAIN_X = linspace(0,length(inversion_at_x_idx_GAIN),length(inversion_at_x_idx_GAIN))*dt;


% TR_ = 30.2730; % measured round trip time;

% middle of the pulse; 
% t0_GAIN = 55.4304;  %  for T1g = 20 (ring cavity/noolap) and p=1.3
% t0_GAIN = 60.5760; % for T1g = 20 (ring cavity/noolap) and p=1.2
[pks,locs ] = findpeaks(E0_sqrt_at_x/E0_sqrt_max,'MinPeakHeight',.6); 
if length(locs)>=2
    t0_GAIN = times_at_GAIN_X(locs(3)); 
else
   t0_GAIN = times_at_GAIN_X(locs(1)); 
end



tau_p_GAIN = FWHM_of_x(x_idx)/1.76;
% tau_p = 2.2625;
T_GAIN = tau_p_GAIN; 
h = 2*tau_p_GAIN/T_GAIN;

t_a = t0_GAIN-TR_/2; t_b = t_a + TR_;
t_l = t0_GAIN-T_GAIN/2; t_r = t0_GAIN+T_GAIN/2; 

dT_GAIN = TR_-T_GAIN;

gamma_1_GAIN = 1/T1_GAIN; 
gamma_2_GAIN = h*R+gamma_1_GAIN;


dt_small  = 30;
t_idx = times_at_GAIN_X >= t_a-dt_small & times_at_GAIN_X <= t_b+dt_small; 
t_vec_new_GAIN = times_at_GAIN_X(t_idx); 
I_vec_new_GAIN = E0_sqrt_at_x(t_idx); 
inv_vec_new_GAIN = inversion_at_x_idx_GAIN(t_idx);


Delta_l_GAIN = d_eq_GAIN*( gamma_1_GAIN/gamma_2_GAIN * (exp(-gamma_2_GAIN*T_GAIN) - 1) - (exp(gamma_1_GAIN*dT_GAIN)-1) )/(exp(-gamma_2_GAIN*T_GAIN)-exp(gamma_1_GAIN*dT_GAIN))
Delta_r_GAIN = (Delta_l_GAIN-d_eq_GAIN* /gamma_2_GAIN)*exp(-gamma_2_GAIN*T_GAIN)+d_eq_GAIN*gamma_1_GAIN/gamma_2_GAIN
Delta_a_GAIN = (Delta_l_GAIN-d_eq_GAIN)*exp(gamma_1_GAIN*dT_GAIN/2)+d_eq_GAIN


boxcar_pulse_GAIN = h*E0_sqrt_max*heaviside(t_vec_new_GAIN - t_l).*heaviside(t_r-t_vec_new_GAIN);




inv_analytical_GAIN = 0*inv_vec_new_GAIN;
d1_idx = t_vec_new_GAIN < t_l;
inv_analytical_GAIN(d1_idx)= (Delta_a_GAIN-d_eq_GAIN)*exp(-gamma_1_GAIN*(t_vec_new_GAIN(d1_idx)-t_a))+d_eq_GAIN; 

d2_idx = t_vec_new_GAIN >= t_l & t_vec_new_GAIN <= t_r;
inv_analytical_GAIN(d2_idx)= (Delta_l_GAIN-d_eq_GAIN*gamma_1_GAIN/gamma_2_GAIN)*exp(-gamma_2_GAIN*(t_vec_new_GAIN(d2_idx)-t_l)) ...
    + d_eq_GAIN*gamma_1_GAIN/gamma_2_GAIN; 

d3_idx = t_vec_new_GAIN > t_r;
inv_analytical_GAIN(d3_idx)= (Delta_r_GAIN-d_eq_GAIN)*exp(-gamma_1_GAIN*(t_vec_new_GAIN(d3_idx)-t_r))+d_eq_GAIN; 

figure;
dTH_vector = d_th*(1+0*t_vec_new_GAIN)
plotyy(t_vec_new_GAIN,[I_vec_new_GAIN;boxcar_pulse_GAIN], ...
    t_vec_new_GAIN,[inv_vec_new_GAIN;inv_analytical_GAIN;dTH_vector]); 


trapz(t_vec_new_GAIN,inv_vec_new_GAIN)
trapz(t_vec_new_GAIN,inv_analytical_GAIN)

%%%% the following is from the analytical expression for the integral! 
int_1 = (Delta_a_GAIN-d_eq_GAIN)/gamma_1_GAIN* ... 
    (1-exp(-gamma_1_GAIN*dT_GAIN/2))+dT_GAIN*d_eq_GAIN/2;

int_2 = 1/gamma_2_GAIN*(Delta_l_GAIN-d_eq_GAIN*gamma_1_GAIN/gamma_2_GAIN)* ...
    (1-exp(-gamma_2_GAIN*T_GAIN))+T_GAIN*d_eq_GAIN*gamma_1_GAIN/gamma_2_GAIN;

int_3 = (Delta_r_GAIN-d_eq_GAIN)/gamma_1_GAIN* ... 
    (1-exp(-gamma_1_GAIN*dT_GAIN/2))+dT_GAIN*d_eq_GAIN/2;

end
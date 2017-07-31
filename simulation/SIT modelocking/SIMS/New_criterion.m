% convert to co-moving frame and plot the field as a function of tau.


% for each tau
close all;

% x_short = x(1:skip:end);
% x_GAIN = x_short(1:size(inversion_data_gain_vs_x,1));
% x_ABS = x_short(size(inversion_data_gain_vs_x,1)+1:end);
% vG = Ltot/TR_;
% 
% time_vec = linspace(-1/2,1/2,size(inversion_data_gain_vs_x,2))*dt*size(inversion_data_gain_vs_x,2);
% 
% x_idx_g = 10;
% x_ = x_GAIN(x_idx_g);
% y_x  = x_ - Ltot/params_gain.L*(x_-params_gain.L);
% [dummy,x_idx_a] = min(abs(x_ABS-y_x));
% 
% 
% tau_vec = linspace(0,2*TR_,1000);
% 
% inversion_tau_g = [];
% inversion_tau_a = [];
% 
% for tau = tau_vec
%     
%     avg_inv = 0;
%     tGAIN_2_take = tau - x_GAIN(x_idx_g)/vG;
%     [dummy,t_idxg] = min(abs(time_vec-tGAIN_2_take));
%     inversion_tau_g = [inversion_tau_g inversion_data_gain_vs_x(x_idx_g,t_idxg)];
%     
%     
%     tABS_2_take  = tau - x_ABS(x_idx_a)/vG;
%     
%     [dummy,t_idxa] = min(abs(time_vec-tABS_2_take));
%     t_idxa;
%     inversion_tau_a = [inversion_tau_a inversion_data_abs_vs_x(x_idx_a,t_idxa)];
% end
% 
% % plotyy(tau_vec,inversion_tau_a,tau_vec,inversion_tau_g);
% 
% Ga = sigma_a*Ncarriers_a*params_abs.L*1e-3;
% Gg = sigma_g*Ncarriers_g*params_gain.L*1e-3;
% powerloss_L = 2*params_gain.linear_loss*100*params_gain.L*1e-3;
% 
% 
% plot(tau_vec,Gg*inversion_tau_g+Ga*inversion_tau_a-powerloss_L);
% hold on
% plot(tau_vec,Gg*inversion_tau_g,tau_vec,Ga*inversion_tau_a);


%%

inversion_tau_g = [];
inversion_tau_a = [];

for tau = tau_vec
    
    avg_inv = 0;
    for j = 1:length(x_GAIN)
        tGAIN_2_take = tau - x_GAIN(j)/vG;
        [dummy,t_idxg] = min(abs(time_vec-tGAIN_2_take));
        avg_inv = avg_inv+inversion_data_gain_vs_x(j,t_idxg);
    end
    inversion_tau_g = [inversion_tau_g avg_inv/length(x_GAIN)];
    
    avg_inv = 0;
    for j = 1:length(x_ABS)
        tABS_2_take  = tau - x_ABS(j)/vG;
        [dummy,t_idxa] = min(abs(time_vec-tABS_2_take));
        avg_inv = avg_inv+inversion_data_abs_vs_x(j,t_idxa);
        
    end
    inversion_tau_a = [inversion_tau_a avg_inv/length(x_ABS)];
end

figure; 
subplot(2,1,1)
plotyy(tau_vec,inversion_tau_a,tau_vec,inversion_tau_g);
subplot(2,1,2)
plot(tau_vec,Gg*inversion_tau_g+Ga*inversion_tau_a-powerloss_L);

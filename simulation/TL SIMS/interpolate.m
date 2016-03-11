
 
clear; clc;  close all;

load('rates'); 
%divide by 10 to normalize to kv/mm; 
vals_fit = vals'/10; 


ac_fit_type = 'poly2';
O_13_fit = fit(vals_fit,O_13_v'/2,ac_fit_type); %in eV

energy_fit_type = 'poly1';

levs = [INJ_1 INJ_2 ULL LLL_1 LLL_2 DEPOP_1 DEPOP_2];

N_lvls = length(levs);
E_fit = cell(length(levs),1); 
for lvl= 1:N_lvls
    E_fit{lvl} = fit(vals_fit,E_v(:,lvl),energy_fit_type); %in eV
end

N_lvls = 7; 

W_fit = cell(N_lvls,N_lvls); 
for lvl_i = 1:N_lvls
    lvl_i
    for lvl_j = 1:N_lvls
        W_fit{lvl_i,lvl_j} = fit(vals_fit,W_v(:,lvl_i,lvl_j)/1E12,'linearinterp'); %in THz
            lvl_j
    end
end

save('fitted_data');

%plot some of the results 
%%
vals_inter = [7:0.2:14]/10;
INJ = 2; ULL = 3; LLL = 4;
dfigure; 
subplot(3,1,1); %plot the energies
hold on;
% plot(vals_inter,E_fit{INJ}(vals_inter),vals_fit,E_v(:,INJ)); 
plot(vals_inter,E_fit{ULL}(vals_inter)-E_fit{LLL}(vals_inter),vals_fit,E_v(:,ULL)-E_v(:,LLL)); 
% plot(vals_inter,E_fit{LLL}(vals_inter),vals_fit,E_v(:,LLL)); 
hold off; 
subplot(3,1,2); %plot the anticrossings
plot(vals_inter,O_13_fit(vals_inter),vals_fit,O_13_v/2); 

subplot(3,1,3); %plot the scattering rates
hold on; 
plot(vals_inter,W_fit{INJ,ULL}(vals_inter),'-r',vals_fit,W_v(:,INJ,ULL)/1E12,'--r'); 
plot(vals_inter,W_fit{ULL,LLL}(vals_inter),'-g',vals_fit,W_v(:,ULL,LLL)/1E12,'--g');
plot(vals_inter,W_fit{ULL,INJ}(vals_inter),'-b',vals_fit,W_v(:,ULL,INJ)/1E12,'--b');
hold off; 


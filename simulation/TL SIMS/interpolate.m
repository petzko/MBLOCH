
 
clear; clc;  

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
        W_fit{lvl_i,lvl_j} = fit(vals_fit,W_v(:,lvl_i,lvl_j)/1E12,'exp2'); %in THz
            lvl_j
    end
end

save('fitted_data');






 
clear; clc;  close all;

load('rates','vals','W_v');
load('HTB');
%divide by 10 to normalize to kv/mm; 
vals_fit_HTB = bias/10; 


ac_fit_type = 'linearinterp';
AC_fit = cell(1,1);
AC_fit{1} = fit(vals_fit_HTB',abs(AC_energies(:,1)),ac_fit_type); %in eV


zUL_fit = fit(vals_fit_HTB',dipoles/10,'linearinterp');

energy_fit_type = 'poly1';

INJ = 1; ULL = 2; LLL = 3; RES = 4; DEPOP= 5; 

levs = [INJ ULL LLL RES DEPOP];

N_lvls = length(levs);
E_fit = cell(N_lvls,1); 

for lvl= 1:N_lvls
    E_fit{lvl} = fit(vals_fit_HTB',Energies(:,lvl),energy_fit_type); %in eV
end
% 
vals_fit_W = vals'/10; 
W_fit = cell(N_lvls,N_lvls); 
for lvl_i = 1:N_lvls
    lvl_i
    for lvl_j = 1:N_lvls
        W_fit{lvl_i,lvl_j} = fit(vals_fit_W,W_v(:,levs(lvl_i),levs(lvl_j))/1E12,'linearinterp '); %in THz
            lvl_j
    end
end

save('fitted_data');

%plot some of the results 
%%
vals_inter = [4:0.2:14]/10;
INJ = 1; ULL = 2; LLL = 3; RES = 4; DEPOP = 5;

dfigure; 
subplot(4,1,1); %plo the energies
hold on;
% plot(vals_inter,E_fit{INJ}(vals_inter),vals_fit_HTB',Energies(:,levs(INJ))); 
% plot(vals_inter,E_fit{ULL}(vals_inter),vals_fit_HTB',Energies(:,levs(ULL))); 
% plot(vals_inter,E_fit{LLL}(vals_inter),vals_fit_HTB',Energies(:,levs(LLL))); 
% plot(vals_inter,E_fit{RES}(vals_inter),vals_fit_HTB',Energies(:,levs(RES))); 
% plot(vals_inter,E_fit{DEPOP}(vals_inter),vals_fit_HTB',Energies(:,levs(DEPOP))); 
plot(vals_fit_HTB',Energies(:,levs(INJ))); 
plot(vals_fit_HTB',Energies(:,levs(ULL))); 
plot(vals_fit_HTB',Energies(:,levs(LLL))); 
plot(vals_fit_HTB',Energies(:,levs(RES))); 
plot(vals_fit_HTB',Energies(:,levs(DEPOP))); 
legend('INJ','ULL','LLL','RES','DEPOP');

% plot(vals_inter,E_fit{INJ}(vals_inter)-E_fit{ULL}(vals_inter)); 
% plot(vals_inter,E_fit{ULL}(vals_inter)-E_fit{LLL}(vals_inter)); 
% plot(vals_inter,E_fit{LLL}(vals_inter)-E_fit{RES}(vals_inter)); 
% plot(vals_inter,E_fit{RES}(vals_inter)-E_fit{DEPOP}(vals_inter)); 
% legend('\Delta E_{INJ,ULL}','\Delta E_{ULL,LLL}','\Delta E_{LLL,RES}','\Delta E_{RES,DEPOP}');

hold off; 

subplot(4,1,2); %plot the anticrossings
plot(vals_inter,AC_fit{1}(vals_inter),vals_fit_HTB,abs(AC_energies(:,1))); 

subplot(4,1,3); %plot the scattering rates
hold on; 
plot(vals_fit_W,W_v(:,levs(INJ),levs(ULL))/1E12,'--r'); 
plot(vals_fit_W,W_v(:,levs(ULL),levs(LLL))/1E12,'--g');
plot(vals_fit_W,W_v(:,levs(ULL),levs(INJ))/1E12,'--b');
plot(vals_fit_W,W_v(:,levs(LLL),levs(DEPOP))/1E12,'-r');
plot(vals_fit_W,W_v(:,levs(RES),levs(DEPOP))/1E12,'-g');
plot(vals_fit_W,W_v(:,levs(LLL),levs(RES))/1E12,'-b');

hold off; 

subplot(4,1,4)
plot(vals_inter,zUL_fit(vals_inter),vals_fit_HTB,dipoles/10)
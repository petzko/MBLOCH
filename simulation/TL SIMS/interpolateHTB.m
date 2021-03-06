
 
clear; clc;  close all;

load('rates','vals','W_v');
load('HTBs');
%divide by 10 to normalize to kv/mm; 
vals_fit_HTB = bias/10; 


ac_fit_type = 'poly2';
AC_fit = cell(2,1);
AC_fit{1} = fit(vals_fit_HTB',abs(AC_energies(:,1)),ac_fit_type); %in eV
AC_fit{2} = fit(vals_fit_HTB',abs(AC_energies(:,2)),ac_fit_type); %in eV

energy_fit_type = 'poly1';
INJ1 = 1; INJ2 = 2; ULL =3 ; LLL1 =4; LLL2 = 5; DEPOP1 = 6; DEPOP2 = 7; 
levs = [INJ1 INJ2 ULL LLL1 LLL2 DEPOP1 DEPOP2];

N_lvls = length(levs); 
E_fit = cell(length(levs),1); 
for lvl= 1:N_lvls
    E_fit{lvl} = fit(vals_fit_HTB',Energies(:,lvl),energy_fit_type); %in eV
end

N_lvls = 7;
vals_fit_W = vals'/10; 
W_fit = cell(N_lvls,N_lvls); 
for lvl_i = 1:N_lvls
    lvl_i
    for lvl_j = 1:N_lvls
        W_fit{lvl_i,lvl_j} = fit(vals_fit_W,W_v(:,lvl_i,lvl_j)/1E12,'linearinterp'); %in THz
            lvl_j
    end
end

save('fitted_data');

%plot some of the results 
%%
vals_inter = [7:0.2:14]/10;
INJ = INJ2; LLL = LLL1;
dfigure; 
subplot(3,1,1); %plot the energies
hold on;
plot(vals_inter,E_fit{INJ}(vals_inter),vals_fit_HTB',Energies(:,INJ)); 
plot(vals_inter,E_fit{ULL}(vals_inter),vals_fit_HTB',Energies(:,ULL)); 
plot(vals_inter,E_fit{LLL}(vals_inter),vals_fit_HTB',Energies(:,LLL)); 
hold off; 
subplot(3,1,2); %plot the anticrossings
hold on
plot(vals_inter,AC_fit{1}(vals_inter),vals_fit_HTB,abs(AC_energies(:,1))); 
plot(vals_inter,AC_fit{2}(vals_inter),vals_fit_HTB,abs(AC_energies(:,2))); 
hold off

subplot(3,1,3); %plot the scattering rates
hold on; 
plot(vals_inter,W_fit{INJ,ULL}(vals_inter),'-r',vals_fit_W,W_v(:,INJ,ULL)/1E12,'--r'); 
plot(vals_inter,W_fit{ULL,LLL}(vals_inter),'-g',vals_fit_W,W_v(:,ULL,LLL)/1E12,'--g');
plot(vals_inter,W_fit{ULL,INJ}(vals_inter),'-b',vals_fit_W,W_v(:,ULL,INJ)/1E12,'--b');
hold off; 


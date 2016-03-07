close all;


n_1 = r11_time;
n_3 = r33_time; 
n_2 = r22_time;
times = linspace(0,length(n_1)*dt,length(n_1));
figure; 
plotyy(times,[n_1;n_3;n_2;pop_time(1,:);pop_time(2,:)],times,abs(E_p).^2);
legend({'n_{INJ2}','n_{ULL}','n_{LLL1}','n_{INJ1}','n_{LL2}','I(t)'})

figure; 
plotyy(times,[n_1-n_3;n_3-n_2],times,abs(E_p).^2);
legend({'n_{INJ2}-n_{ULL}','n_{ULL}-n_{LLL}','I(t)'})

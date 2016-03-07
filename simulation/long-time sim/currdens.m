close all;
J = [];
Lp =54.8*1E-9; %in meters!   
rt_start = 600; rt_end = 995;
iter_per_rt = round(T_R/dt/iterperrecord);

rates = r11_time(rt_start*iter_per_rt+1:rt_end*iter_per_rt)*W(INJ,DEPOP) + r22_time(rt_start*iter_per_rt+1:rt_end*iter_per_rt)*W(LLL,DEPOP) +...
    r33_time(rt_start*iter_per_rt+1:rt_end*iter_per_rt)*W(ULL,DEPOP);
for p = 1:N_rest
    p_glob_idx = global_idx_rest(p);
    rates = rates + pop_time(p,(rt_start*iter_per_rt+1:rt_end*iter_per_rt))*W(p_glob_idx,DEPOP);
end

%in A/cm^2
Jm2 = Constants('q0')*Lp*Ncarriers/trace_rho*1E12*rates;  J = Jm2/(1E4); tm =  linspace(0,tEnd,length(rates));
subplot(2,1,1); plot(tm,J/1E3); xlabel('time [ps]'); ylabel('curr densiy kA/cm^2'); xlim([tm(1),tm(end)]);
apodize = hanning(length(J))';
J_w = fft(J.*apodize); J_w = J_w(1:end/2-1); f_J = linspace(0,1/2,length(J_w))/dt;
subplot(2,1,2); plot(f_J, abs(J_w)); xlabel('Freq. [THz]'); ylabel('Microwave signal (a.u.)'); xlim([0.1*f_R,10*f_R]);


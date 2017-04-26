clc;clear;close all;

%*************************************************************************%
%****************** Set up for this script *******************************%
% Give folder name of sim files here !!!
simFolder = '4.25';
RT1 = 391; RT2 = 395;    % Plot pulses in time domain
RT3 = 200;               % FFT the last 200 round trips
select_mod = 1;          % Select legend for modA(1) or modF(2)
%*************************************************************************%
%*************************************************************************%



%********************** Load parameters **********************************%
currentFolder = pwd;
FolderPath = fullfile(currentFolder,'/',simFolder);
userpath(FolderPath); % add sim folder into user path

simFiles = dir ([simFolder,'/*.mat']); 
N = length(simFiles); % amount of *.mat files
modA = zeros(N,1);
modF = zeros(N,1);
bias = zeros(N,1);
formatSpec1 = 'modA = %.1f V';
formatSpec2 = 'f_{RF} = %.1f GHz';


for k = 1:N
    
    sim{k} = matfile(simFiles(k).name);
    sim_params{k} = sim{k}.sim_params;
    modA(k) = sim_params{k}.modA;    % modulation amplitude
    modF(k) = sim_params{k}.modF;    % modulation frequency unit: THz
    bias(k) = sim_params{k}.bias;    % bias unit: kV/cm
    strA{k} = sprintf(formatSpec1,modA(k));
    strF{k} = sprintf(formatSpec2,modF(k)*1000);
end

dt = sim{1}.dt;                    % time step unit: ps
T_R = sim{1}.T_R;                  % round trip time unit: ps
record_time = sim{1}.record_time;  % time axis
length_QCL = sim_params{1}.length;  % QCL length unit: mm
% R, L, G, C
rlgc = sim{1}.rlgc;      
rlgc.C = rlgc.C*1e-12;      % Distributed capacitance unit: F/mm
rlgc.L = rlgc.L*1e-12;      % Distributed inductance  unit: H/mm
rlgc.G = 3e-2;              % Distributed conductance unit: S/mm
rlgc.R = 45*sqrt(modF)*length_QCL; % QCL resistance unit: Ohm --from paper W. Maineult: Microwave modulation of terahertz quantum cascade lasers: a transmission-line approach
Res_TL = zeros(N,1);
for k = 1:N
    Res_TL(k) = real(sqrt((rlgc.R(k)+1i*2*pi*modF(k)*1e12*rlgc.L)/(rlgc.G+1i*2*pi*modF(k)*1e12*rlgc.C)));
end

f_R = 1/T_R;                       % round trip frequency unit: THz
rows = length(record_time);
t1 = floor(RT1*T_R/dt);
t2 = floor(RT2*T_R/dt);
%********************  Display of legends ********************************%

if select_mod == 1
    str = strA;
else str = strF;
end
%***************** Load datas from each sim file *************************%
record_U = zeros(rows,N);
record_V = zeros(rows,N);
record_v_TL = zeros(rows,N);
record_i_TL = zeros(rows,N);
record_J_TL = zeros(rows,N);
ac_I = zeros(N,1);      % fraction of ac current component
ac_U = zeros(N,1);      % fraction of ac voltage component

for k = 1:N

    record_U(:,k) = sim{k}.record_U;
    record_V(:,k) = sim{k}.record_V;
    record_v_TL(:,k) = sim{k}.record_v_TL;
    record_i_TL(:,k) = sim{k}.record_i_TL;
    record_J_TL(:,k) = sim{k}.record_J_TL;
end

userpath('clear');  % clear user path



%*************************************************************************%
%************************** Pulses in 5 round trips **********************%
x = dt*linspace(t1,t2,t2-t1+1)/T_R;
Ez_5RT = zeros(length(x),N);
for k = 1:N
    Ez_5RT(:,k) = abs(record_U(t1:t2,k));
    Ez_5RT(:,k) = Ez_5RT(:,k)/max(Ez_5RT(:,k));
end
%*************************************************************************%
%************************* Freuquency spectra*****************************%
index = floor(RT3*T_R/dt);
Ez_abs = zeros(index,N);
Ez_shift = zeros(index,N);
taxis = dt*linspace(rows-index+1,rows,index);
taxis = taxis';

for k = 1:N
    Ez_abs(:,k) = abs(record_U(end-index+1:end,k));
    Ez_shift(:,k) = Ez_abs(:,k).*exp(1i*2*pi*3.8*taxis);
    [ FFT,f1 ] = mydft2(Ez_shift(:,k),dt);
    if k == 1
       Ez3 = zeros(length(FFT),N);
    end
    Ez3(:,k) = abs(FFT(1:length(FFT)));
    Ez3(:,k) = Ez3(:,k)/max(Ez3(:,k));
end

f1 = f1'/1E12;  f1 = f1(1:length(f1));  % Frequency unit: THz
%*************************************************************************%
%************************* Beatnote **************************************%
for k = 1:N
    [ FFT,f2 ] = mydft2(Ez_abs(:,k),dt);
    if k == 1
       Ez4 = zeros(length(FFT)/2,N);
    end
    Ez4(:,k) = abs(FFT(1:length(FFT)/2));
    Ez4(:,k) = Ez4(:,k)/max(Ez4(:,k));
    

    [ FFT,~ ] = mydft2(record_v_TL(end-index+1:end,k),dt);
    if k == 1
       v1 = zeros(length(FFT)/2,N);
    end
    v1(:,k) = abs(FFT(1:length(FFT)/2));
    v1(:,k) = v1(:,k)/max(v1(100:end,k));
    
end

f2 = f2'/1E9;  f2 = f2(1:length(f2)/2);  % Frequency unit: THz
%*************************************************************************%
%************************* RF power calculation **************************%
for k = 1:N
    % ac current fraction
    [FFT, f3] = mydft2(record_i_TL(end/2:end,k),dt);
    if k == 1
       P_i = zeros(length(FFT),N);
       P_v = P_i;
       indices = find(abs(f3-1.3e10)<6e9);
    end
    P_i(:,k) = abs(FFT(1:length(FFT)));
    ac_I(k) = max(P_i(indices,k))/max(P_i(:,k));
    % ac votage fraction
    [FFT, ~] = mydft2(record_v_TL(end/2:end,k),dt);
    P_v(:,k) = abs(FFT(1:length(FFT)));
    ac_U(k) = max(P_v(indices,k))/max(P_v(:,k));

end
ac_P = 10*log10(bias.*ac_U.*ac_I*record_i_TL(end,1)*1e3);
ac_P2 = 10*log10(Res_TL.*(ac_I*record_i_TL(end,1)).^2*1e3);
display('RF power (dBm)');
display(ac_P);
display(ac_P2);
%*************************************************************************%
%************************* Pulse width calculation ***********************%
t3 = floor((RT1+1)*T_R/dt); % One round trip, single pulse
FWHM = zeros(N,1);
Ez_t = zeros(t3-t1+1,N);

for k=1:N
Ez_t(:,k) = record_U(t1:t3,k);
Ez_t(:,k) = abs(Ez_t(:,k));
Ez_t(:,k) = Ez_t(:,k)/max(Ez_t(:,k));

t_pulse = Ez_t(Ez_t(:,k)>=0.5);
FWHM(k) = length(t_pulse)*dt;
end
x2 = dt*linspace(t1,t3,t3-t1+1);
display('Pulse width (ps)');
display(FWHM);

% figure(4)
% plot(x,Ez_t)
% xlabel('Time (ps)')
% ylabel('Normalized intensity (a.u.)')
% text(4655.4,0.5,'\leftarrow 11.6ps','Color','red','FontSize',12)
% text(4640.8,0.5,'\rightarrow','Color','red','FontSize',12)


%*************************************************************************%
%************************ Plot figures ***********************************%
%*************************************************************************%

% Figure 1: Pulses in time domain with 5 round trips
figure(1)
for k = 1:N
    subplot(N,1,k)
    plot(x',Ez_5RT(:,k))
    ylim([0,1.2])
    xlim([RT1,RT2])
    if k == 1
        if select_mod == 2
            title_name = sprintf('Pulses in 5 round trips (t_{rt} = 72.0498 ps, modA = %.1f V, bias = %.1f kV/cm)',modA(1), bias(1));
        else title_name = sprintf('Pulses in 5 round trips (t_{rt} = 72.0498 ps, f_{RF} = %.1f GHz)',modF(1)*1000);        
        end
        title(title_name)
    end
    ylabel('normalized |E_z|')
    if k == N
        xlabel('Round trip')
    end

    text(RT2-0.8,0.8,str(k),'Color','red','FontSize',10)

end

% Figure 2: Frequency spectra
figure(2)
for k = 1:N
    subplot(N,1,k)
    plot(f1,Ez3(:,k))
    xlim([3.5,4.1])
    ylim([0,1.2])
    if k == 1
        if select_mod == 2
            title_name = sprintf('Frequency spectra (modA = %.1f V, bias = %.1f kV/cm)',modA(1),bias(1));
        else title_name = sprintf('Frequency spectra (f_{RF} = %.1f GHz)',modF(1)*1000);        
        end
        title(title_name)
    end
    ylabel('Intencity a.u.')
    if k == N
        xlabel('Frequency /THz')
    end
    text(4,0.8,str(k),'Color','red','FontSize',10)

end

% Figure 3: Beatnote
figure(3)
for k = 1:N
    subplot(N,1,k)
    yyaxis left
    semilogy(f2,Ez4(:,k))
    xlim([2,18])
    ylim([1e-8,10])
    if k == 1
        if select_mod == 2
            title_name = sprintf('Beatnote (modA = %.1f V, bias = %.1f kV/cm)',modA(1),bias(1));
        else title_name = sprintf('Beatnote (f_{RF} = %.1f GHz)',modF(1)*1000);
        end
        title(title_name)
    end
    ylabel('Intencity /dB')
    if k == N
        xlabel('Frequency /GHz')
    end
    text(15.5,0.001,str(k),'Color','red','FontSize',10)
    
    
    yyaxis right
    semilogy(f2,v1(:,k))
    xlim([2,18])
    ylim([1e-8,10])
    legend('E_z','v\_TL','Location','northwest')
    
end
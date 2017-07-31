close all; 
Length = 5*1E-3; % 5mm
Width = 40*1E-6; % 20um 
Height = 10e-6; 
envelope = record_U_a;
n = 3.6;
tch=1e-12;
zUL = sim_settings.zULa
dt = dat.dt;
normalization_factor = 1E12*Constants('hbar')/(Constants('q0')*zUL*1E-9);
I_t = Constants('eps0')*n*Constants('c')/8*(abs((envelope)*normalization_factor).^2); % in watts
power = I_t*Width*Height*1e3;
time =linspace(0, dt*length(I_t),length(I_t))*tch;

plot(time,power);
xlabel('time [s]');
ylabel('intensity [mW]');





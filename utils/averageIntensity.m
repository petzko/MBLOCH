Length = 5*1E-3; % 5mm
Width = 20*1E-6; % 20um 
Height = 183*54.8E-9; 
a_m = .2; %mirror losses for a 2mm cavity 
R = .80

normalization_factor = 1E12*Constants('hbar')/(Constants('q0')*zUL*1E-9);
I_t = Constants('eps0')*n*Constants('c')/8*(abs((E_p+E_m)*normalization_factor).^2);
inPower = I_t*Width*Height*1E3; %miliwatts

time =linspace(0, dt*length(I_t),length(I_t))*tch;

avgInPower = mean(inPower)
maxInPower = max(inPower)
outPower = (1-R)/R*inPower; %in miliwats
plot(time,outPower);
xlabel('time [s]');
ylabel('output power [mW]');

avgOutPower = mean(outPower)
maxOutPower = max(outPower)


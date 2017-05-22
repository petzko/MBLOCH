%Due to symmetry of the transmission line, the scattering matrix
% can be simplified with S11 = S22 and S12 = S21. 
Comsol=[ 3.0000	-0.85641+0.11367i	0.067896-0.069586i ];  %from Comsol

S11 = Comsol(2);
S21 = Comsol(3);
S22=S11;
S12=S21;
s_params = [S11,S12; S21,S22];

length = 2.65e-3;       % cavity length 2.65mm
freq = Comsol(1)*1e9;   % mode frequency
z0 = 20.977;

rlgc_params = s2rlgc(s_params,length,freq,z0);
sprintf('Frequency @ %0.2f GHz', Comsol(1))
display(rlgc_params);	
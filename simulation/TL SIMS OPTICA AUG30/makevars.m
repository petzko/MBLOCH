load('rates'); 

W_THz = W_v/1E12; 
BIAS = vals; 


load('HTB')
BIAS02 = bias(16:end); 

ENERGIES = Energies(16:end,:); %in eV
ANTICROSSINGS = AC_energies(16:end,1); %in eV 
DIPOLES = dipoles(16:end,:)/10; %in NM



save('rawdata','W_THz','BIAS','ENERGIES','ANTICROSSINGS','DIPOLES','INJ','ULL','LLL','RES','DEPOP');

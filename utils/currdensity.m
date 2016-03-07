files = dir(pwd);
currCtr  =1;
I = [];
J = [];
Length = 5*1E-3;% 5mm
Width = 20*1E-6;% 20 um 
Lp =54.8*1E-9; %in meters!   
        
for i = 1:length(files)
    
    f  = files(i);
    if (strfind(f.name,'.mat'))
        f.name
        load(f.name);
        Jm2 = Constants('q0')*Lp*Ncarriers/trace_rho*1E12*(mean(r110)*W(INJ,DEPOP) + mean(r220)*W(LLL,DEPOP) + mean(r330)*W(ULL,DEPOP)); 
        J(currCtr) = Jm2/(1E4); %in A/cm^2
        I(currCtr) = Jm2*Length*Width; %in A
        currCtr =currCtr +1;
    end
    
end
J
I
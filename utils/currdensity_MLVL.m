files = dir(pwd);
currCtr  =1;
I = [];
J = [];
Length = 5*1E-3;% 5mm
Width = 20*1E-6;% 20 um 
Lp =54.8*1E-9; %in meters!   
        pattern = 'MULTILEVELS_.+\.mat'
for i = 1:length(files)
    
    f  = files(i);
    if (~isempty(regexp(f.name,pattern)))
        f.name
        load(f.name);
        rates = mean(r110)*W(INJ,DEPOP) + mean(r220)*W(LLL,DEPOP) + mean(r330)*W(ULL,DEPOP);
        for p = 1:N_rest
            p_glob_idx = idx_rest(p);
            rates = rates + mean(populations(:,p))*W(p_glob_idx,DEPOP);
        end
        
        Jm2 = Constants('q0')*Lp*Ncarriers/trace_rho*1E12*rates; 
        J(currCtr) = Jm2/(1E4); %in A/cm^2
        I(currCtr) = Jm2*Length*Width; %in A
        currCtr =currCtr +1;
    end
    
end
J
I
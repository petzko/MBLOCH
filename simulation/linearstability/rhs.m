function [res] = rhs(in,params)
    % idx =[ 1    2   3    4    5    6   7    8     9   10    11    12    13    14    15    16          17]
    % in = [f_p,f_m,r110,r330,r220,r11p,r33p,r22p,r130,r13_p,r13_m,n12_p,n12_m,n32_p,n32_m,pop_INJ1,pop_DEPOP1]
    INJ = params.INJ; ULL = params.ULL; LLL = params.LLL; DEPOP= params.DEPOP; INJ2 = params.INJ2; DEPOP2 = params.DEPOP2; W= params.W; 
    O13 = params.O13; G = params.G; dE13 = params.dE13; dE32 = params.dE32; dE12 = params.dE12; 
    %%%%%%%%%%%%%%%%%%%%%%%%
    losses = params.losses;
    
    %%%%%% fields 
    res(1) = -1i*c*in(14)-c*losses*in(1);
    res(2) = -1i*c*in(15)-c*losses*in(2);
    
    %%%%%% BLOCH PART %%%%%%
    %%%% POPULATIONS
    %%% r110 , 3 
    res(3) = 1i*params.O13.*(in(9)-conj(in(9))) +( params.W(ULL,INJ) + params.W(ULL,DEPOP) )*in(4)+ (W(LLL,INJ)+W(LLL,DEPOP))*in(5)  - G(INJ)*in(3);
    res(3) = res(3)+ (W(INJ2,INJ)+W(INJ2,DEPOP))*in(16)+(W(DEPOP2,INJ)+W(DEPOP2,DEPOP))*in(17);
    
    %%% r330 , 4
    lmInteraction = conj(in(1)).*in(14) + conj(in(2)).*in(15);
    res(4) = 1i*O13.*(conj(in(9)) - in(9)) +1i/2.*(lmInteraction-conj(lmInteraction)) + in(3)*W(INJ,ULL) + in(5)*W(LLL,ULL) - G(ULL)*in(4);
    res(4) = res(4)+W(INJ2,ULL)*in(16)+W(DEPOP2,ULL)*in(17);
    
    %%% r220 , 5
    res(5) = -1i/2.*(lmInteraction-conj(lmInteraction)) + in(3)*W(INJ,LLL) + in(4)*W(ULL,LLL) - G(LLL)*in(5);
    res(5) = res(5)+W(INJ2,LLL)*in(16)+W(DEPOP2,LLL)*in(17);
    
    %%% rINJ2 and rDEPOP2 16,17
    res(16) = W(INJ,INJ2)*in(3)+W(ULL,INJ2)*in(4)+W(LLL,INJ2)*in(5) + W(DEPOP2,INJ2)*in(17) - G(INJ2)*in(16);
    res(17) = W(INJ,DEPOP2)*in(3)+W(ULL,DEPOP2)*in(4)+W(LLL,DEPOP2)*in(5) + W(INJ2,DEPOP2)*in(16) - G(DEPOP2)*in(17);
    
    %%%% COHERENCES
    %%% r130 - 9
    res(9) = dE13*in(9) + 1i*O13*(in(3)-in(4)) + 1i/2*(conj(in(1)).*in(12) + conj(in(2)).*in(13));
    
    %%% n32p - 14
    res(14) = dE32*in(14) + 1i/2*(in(1).*(in(4)-in(5)) + in(2).*(in(7)-in(8))) - 1i*O13*in(12);
    
    %%% n32m = 15; 
    res(15) = dE32*in(15) + 1i/2*(in(2).*(in(4)-in(5)) + in(1).*conj(in(7)-in(8))) - 1i*O13*in(13);
    
    %%% n12p = 12; 
    res(12) = dE12*in(12) +1i/2*(in(1).*in(9) + in(2).*in(10)) - 1i*O13*in(14);

    %%% n12m = 13; 
    res(13) = dE12*in(13) + 1i/2*(in(2).*in(9) +in(1).*in(11)) - 1i*O13*in(15);
    
  
    % idx =[ 1    2   3    4    5    6   7    8     9   10    11    12    13    14    15    16          17]
    % in = [f_p,f_m,r110,r330,r220,r11p,r33p,r22p,r130,r13_p,r13_m,n12_p,n12_m,n32_p,n32_m,pop_INJ1,pop_DEPOP1]
        
    %%% r11+ 6
    res(6) = 1i*O13.*(in(10)-conj(in(11))) + (W(ULL,INJ)+W(ULL,DEPOP))*in(7)+ (W(LLL,INJ)+W(LLL,DEPOP))*in(8) - (G(INJ)+diffusion)*in(6);
        
    %%% r33+ 7
    res(7) = 1i*O13.*(conj(in(11))-in(10))+1i/2*(conj(in(2)).*(in(14)) -(in(1)).*conj(in(15))) +  W(INJ,ULL)*in(6) + W(LLL,ULL)*in(8) - (G(ULL)+diffusion)*in(7);
   
    %%% r22+ 8
    res(8) = -1i/2*(conj(in(2)).*(in(14)) -(in(1)).*conj(in(15))) +  W(INJ,LLL)*in(6) + W(ULL,LLL)*in(7) - (G(LLL)+diffusion)*in(8);
        
    %%% r13+ 10
    res(10) = (dE13-diffusion)*in(10) + 1i*O13*(in(6)-in(7)) +1i/2*(conj(in(2)).*in(12));
        
    %%% r13- 11
    res(11) = (dE13-diffusion)*in(11) + 1i*O13*conj(in(6)-in(7)) + 1i/2*(conj(in(1)).*in(13));
    

end
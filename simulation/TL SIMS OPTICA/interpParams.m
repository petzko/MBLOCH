function dat = interpParams(settings,dat,W_fit,E_fit,AC_fit,zUL_fit)

  %%% begin setting up TL params!
        INJ = dat.INJ; ULL = dat.ULL; LLL = dat.LLL; RES = dat.RES; DEPOP = dat.DEPOP; 
        
        dat.W(:,INJ,ULL) = W_fit{INJ,ULL}(dat.v_TL);
        dat.W(:,INJ,LLL) = W_fit{INJ,LLL}(dat.v_TL);
        dat.W(:,INJ,RES) = W_fit{INJ,RES}(dat.v_TL);
       
        dat.W(:,ULL,INJ) = W_fit{ULL,INJ}(dat.v_TL);
        dat.W(:,ULL,LLL) = W_fit{ULL,LLL}(dat.v_TL);
        dat.W(:,ULL,RES) = W_fit{ULL,RES}(dat.v_TL);
        dat.W(:,ULL,DEPOP) = W_fit{ULL,DEPOP}(dat.v_TL);
%         
%         dat.W(:,LLL,INJ) = W_fit{LLL,INJ}(dat.v_TL);
%         dat.W(:,LLL,ULL) = W_fit{LLL,ULL}(dat.v_TL);
%         dat.W(:,LLL,RES) = W_fit{LLL,RES}(dat.v_TL);
%         dat.W(:,LLL,DEPOP) = W_fit{LLL,DEPOP}(dat.v_TL);
% 
%         dat.W(:,RES,INJ) = W_fit{RES,INJ}(dat.v_TL);
%         dat.W(:,RES,ULL) = W_fit{RES,ULL}(dat.v_TL);
%         dat.W(:,RES,LLL) = W_fit{RES,LLL}(dat.v_TL);
%         dat.W(:,RES,DEPOP) = W_fit{RES,DEPOP}(dat.v_TL);
%         
%         dat.G(:,INJ) = dat.W(:,INJ,ULL)+ dat.W(:,INJ,LLL)+dat.W(:,INJ,RES);
%         dat.G(:,ULL) = dat.W(:,ULL,INJ)+dat.W(:,ULL,LLL)+dat.W(:,ULL,DEPOP) + dat.W(:,ULL,RES);
%         dat.G(:,LLL) =  dat.W(:,LLL,INJ)+dat.W(:,LLL,ULL)+dat.W(:,LLL,DEPOP) + dat.W(:,LLL,RES) ;
%         dat.G(:,RES) =  dat.W(:,RES,INJ)+dat.W(:,RES,ULL)+dat.W(:,RES,LLL) +dat.W(:,RES,DEPOP) ;

        %%% begin setting up TL params!
        dat.G(:,INJ) = 0*dat.G(:,INJ);
        for el = setdiff(1:dat.NLVLS,[dat.INJ,dat.DEPOP])
            dat.W(:,INJ,el) = W_fit{INJ,el}(dat.v_TL);
            dat.G(:,INJ) = dat.G(:,INJ)+dat.W(:,INJ,el);
        end

        for lvl = setdiff(1:dat.NLVLS,[dat.INJ,dat.DEPOP])
            dat.G(:,lvl) = 0* dat.G(:,lvl);
            for el = 1:dat.NLVLS
                if el ~=lvl
                    dat.W(:,lvl,el) = W_fit{lvl,el}(dat.v_TL);
                    dat.G(:,lvl) = dat.G(:,lvl) +  dat.W(:,lvl,el);
                end
            end
        end

        
        % % % % Added Pure Dephasing
        if(settings.deph>0)
            dat.gamma_13 = 0.5*(dat.G(:,dat.INJ)+dat.G(:,dat.ULL))+1/settings.Tdeph_1; %% dephsing of the resonant tunneling transition
            dat.gamma_32 = 0.5*(dat.G(:,dat.ULL)+dat.G(:,dat.LLL))+1/settings.Tdeph_2; % dephasing of the optical transision...
            dat.gamma_12 = 0.5*(dat.G(:,dat.INJ)+dat.G(:,dat.LLL))+1/settings.Tdeph_3; % dephasing of the latest transition
        else
            dat.gamma_13 = 0.5*(dat.G(:,dat.INJ)+dat.G(:,dat.ULL)); %% dephsing of the resonant tunneling transition
            dat.gamma_32 = 0.5*(dat.G(:,dat.ULL)+dat.G(:,dat.LLL)); % dephasing of the optical transision...
            dat.gamma_12 = 0.5*(dat.G(:,dat.INJ)+dat.G(:,dat.LLL)); % dephasing of the latest transition
        end
        
         %%%%% correctly setup the energies, the anticrossings and the scattering rates...
        %obtain the energies of the core levels for the simulation
        E1 = E_fit{dat.INJ}(dat.v_TL)/dat.hbar;
        E3 = E_fit{dat.ULL}(dat.v_TL)/dat.hbar;
        E2 = E_fit{dat.LLL}(dat.v_TL)/dat.hbar; % in rad/ps
    

        %%%% energies: E1 > E3 > E2 Optical transition is 3->2 TUNNELING transition
        %%%% is 1->3
        %%%% energies: E1 > E3 > E2 Optical transition is 3->2 TUNNELING transition is 1->3
        dat.O13 = AC_fit{1}(dat.v_TL)/dat.hbar; % in rad/ps

        dat.E13 = E1-E3; %rad/ps; 1->3 traisition freq
        dat.E12 = E1-E2; %rad/ps; 1->2 transition freq
        dat.E32 = E3-E2; %rad/ps  3->2 transition freq (optical transition)

        dat.dE13 = -1i*dat.E13 - dat.gamma_13; %
        dat.dE32 = +1i*(dat.E0 - dat.E32) - dat.gamma_32; %
        dat.dE12 = +1i*(dat.E0 - dat.E12)- dat.gamma_12; %
        
        %the varying dipole ratio. 
        dat.dipR = zUL_fit(dat.v_TL)/dat.zUL;
        
        
      
end
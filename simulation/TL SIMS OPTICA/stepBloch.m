function dat = stepBloch(settings,dat)
    
    W = dat.W; INJ = dat.INJ; ULL = dat.ULL; LLL = dat.LLL; RES = dat.RES; DEPOP = dat.DEPOP; G = dat.G;
    %time resolved dipole moment ratio.
    
    %%%% POPULATIONS
    r110_t = 1i*dat.O13.*(dat.r130-conj(dat.r130)) +(W(:,ULL,INJ) + W(:,ULL,DEPOP)).*dat.r330 + (W(:,LLL,INJ)+W(:,LLL,DEPOP)).*dat.r220 + (W(:,RES,INJ)+W(:,RES,DEPOP)).*dat.rRES- G(:,INJ).*dat.r110;
    dat.r110_solver.make_step(r110_t,dat.dt);
    
    lmInteraction = dat.dipR.*(conj(dat.U).*dat.n32p + conj(dat.V).*dat.n32m);
    r330_t = 1i*dat.O13.*(conj(dat.r130) - dat.r130) +1i/2*(lmInteraction-conj(lmInteraction)) +  dat.r110.*W(:,INJ,ULL) + dat.r220.*W(:,LLL,ULL) + dat.rRES.*W(:,RES,ULL)- G(:,ULL).*dat.r330;
    dat.r330_solver.make_step(r330_t,dat.dt);
    
    r220_t = -1i/2*(lmInteraction-conj(lmInteraction)) + dat.r110.*W(:,INJ,LLL)+ dat.r330.*W(:,ULL,LLL) +dat.rRES.*W(:,RES,LLL)- G(:,LLL).*dat.r220;
    dat.r220_solver.make_step(r220_t,dat.dt);
    
    rRES_t =  dat.r110.*W(:,INJ,RES)+ dat.r330.*W(:,ULL,RES)+ dat.r220.*W(:,LLL,RES) - G(:,RES).*dat.rRES;
    dat.rRES_solver.make_step(rRES_t,dat.dt);
    
    %%%% COHERENCES
    r130_t = dat.dE13.*dat.r130 + 1i*dat.O13.*(dat.r110-dat.r330) +1i/2*dat.dipR.*(conj(dat.U).*dat.n12p + conj(dat.V).*dat.n12m);
    dat.r130_solver.make_step(r130_t,dat.dt);
    
    dat.n32p_t = dat.dE32.*dat.n32p + 1i/2*dat.dipR.*(dat.U.*(dat.r330-dat.r220) + dat.V.*(dat.r33p-dat.r22p)) - 1i*dat.O13.*dat.n12p;
    dat.n32p_solver.make_step(dat.n32p_t,dat.dt);
    
    dat.n32m_t = dat.dE32.*dat.n32m + 1i/2*dat.dipR.*(dat.V.*(dat.r330-dat.r220) + dat.U.*conj(dat.r33p-dat.r22p)) - 1i*dat.O13.*dat.n12m;
    dat.n32m_solver.make_step(dat.n32m_t,dat.dt);
    
    n12p_t = dat.dE12.*dat.n12p +1i/2*dat.dipR.*(dat.U.*dat.r130 + dat.V.*dat.r13p) - 1i*dat.O13.*dat.n32p;
    dat.n12p_solver.make_step(n12p_t,dat.dt);
    
    n12m_t = dat.dE12.*dat.n12m + 1i/2*dat.dipR.*(dat.V.*dat.r130 +dat.U.*dat.r13m) - 1i*dat.O13.*dat.n32m;
    dat.n12m_solver.make_step(n12m_t,dat.dt);
    
    if(settings.shb > 0 )
        
        %%% r11+
        r11p_t = 1i*dat.O13.*(dat.r13p-conj(dat.r13m)) + (W(:,ULL,INJ)+W(:,ULL,DEPOP)).*dat.r33p+ (W(:,LLL,INJ)+W(:,LLL,DEPOP)).*dat.r22p - (G(:,INJ)+dat.diffusion).*dat.r11p;
        dat.r11p_solver.make_step(r11p_t,dat.dt);
        
        %%% r33+
        r33p_t = 1i*dat.O13.*(conj(dat.r13m)-dat.r13p)+1i/2*dat.dipR.*(conj(dat.V).*(dat.n32p) -(dat.U).*conj(dat.n32m)) +  W(:,INJ,ULL).*dat.r11p + W(:,LLL,ULL).*dat.r22p - (G(:,ULL)+dat.diffusion).*dat.r33p;
        dat.r33p_solver.make_step(r33p_t,dat.dt);
        
        %%% r22+
        r22p_t = -1i/2*dat.dipR.*(conj(dat.V).*(dat.n32p) -(dat.U).*conj(dat.n32m)) + W(:,INJ,LLL).*dat.r11p + W(:,ULL,LLL).*dat.r33p - (G(:,LLL)+dat.diffusion).*dat.r22p;
        dat.r22p_solver.make_step(r22p_t,dat.dt);
        
        %%% r13+
        r13p_t = (dat.dE13-dat.diffusion).*dat.r13p + 1i*dat.O13.*(dat.r11p-dat.r33p) +1i/2*dat.dipR.*(conj(dat.V).*dat.n12p);
        dat.r13p_solver.make_step(r13p_t,dat.dt);
        
        %%% r13-
        r13m_t = (dat.dE13-dat.diffusion).*dat.r13m + 1i*dat.O13.*conj(dat.r11p-dat.r33p) + 1i/2*dat.dipR.*(conj(dat.U).*dat.n12m);
        dat.r13m_solver.make_step(r13m_t,dat.dt);
    end
    
     dat.rates =   dat.r110.*dat.W(:,dat.INJ,dat.DEPOP) + dat.r220.*dat.W(:,dat.LLL,dat.DEPOP) + ...
    dat.r330.*dat.W(:,dat.ULL,dat.DEPOP)+dat.rRES.*dat.W(:,dat.RES,dat.DEPOP);
    dat.rates_t =  r110_t.*dat.W(:,dat.INJ,dat.DEPOP) + r220_t.*dat.W(:,dat.LLL,dat.DEPOP) + ...
    r330_t.*dat.W(:,dat.ULL,dat.DEPOP)+rRES_t.*dat.W(:,dat.RES,dat.DEPOP);
    

end
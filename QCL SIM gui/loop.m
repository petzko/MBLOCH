

while(t< tEnd)
   
    %%plot some of the results if neeed ariseth :D
    
    %%plot some of the results if neeed ariseth :D
    if(mod(iter_ctr,plotCtr) == 0)
        clc;
        display(['iteration Nr. = ' num2str(iter_ctr) ' @ RT = ' num2str(t/T_R)])
        display(['outcoupl: ' num2str(out_coupling)])
        display(['sat abs coeff: ' num2str(min(l))])
        intensity = pField.*conj(pField) + mField.*conj(mField);
        maxInt = max(intensity);
        display(['max Intensity: ' num2str(maxInt) ]);
% %         ax = plotyy(x,intensity,x,[r110,(r220),(r330)]);
% %         set(ax(1),'Xlim',[-0.2 Ltot+0.2]);
% %         set(ax(2),'Xlim',[-0.1 Ltot+0.2]);
% %         title(['IME + RNFD  @ t = ' num2str(t)]);
        getframe;
    end
    
    %%%% obtain the field, field intensity and the total population at
    %%%% position "idx" ...
    if(t >= tEnd - pulse_len && mod(iter_ctr,recordIter) == 0)
        
        intensity = pField(idx).*conj(pField(idx)) + mField(idx).*conj(mField(idx));
        I_t(ctr) = intensity;
        E_p(ctr)= pField(idx);
        E_m(ctr)= mField(idx); 
        r11_t(ctr)= r110(idx);
        r33_t(ctr)= r330(idx);
        r22_t(ctr) = r220(idx);
        pop_sum(ctr) = r110(idx)+r220(idx)+r330(idx);
        ctr = ctr+1;   
        
    end
    
    %%%%%% BLOCH PART %%%%%%
    %%%% POPULATIONS
    r110_t = 1i*O13.*(conj(r310)-r310) + r330*W31 + r220*W21 - G1*r110;
    r110_solver.make_step(r110_t,dt);
    
    lmInteraction = pField.*conj(n23p) + mField.*conj(n23m);
    r330_t = 1i*O13.*(r310 - conj(r310)) +1i/2.*(lmInteraction-conj(lmInteraction)) + r110*W13 + r220*W23 - G3*r330;
    r330_solver.make_step(r330_t,dt);
    
    r220_t = -1i/2.*(lmInteraction-conj(lmInteraction)) + r110*W12 + r330*W32 - G2*r220;
    r220_solver.make_step(r220_t,dt);
   
    %%%% COHERENCES
    r310_t = dE31*r310 - 1i*O13*(r110-r330) -1i/2*(conj(pField).*n21p + conj(mField).*n21m);
    r310_solver.make_step(r310_t,dt);
   
    n23p_t = dE23*n23p - 1i/2*(pField.*(r330-r220) + mField.*conj(r33p-r22p)) + 1i*O13*n21p; 
    n23p_solver.make_step(n23p_t,dt); 
    
    n23m_t = dE23*n23m - 1i/2*(mField.*(r330-r220) + pField.*(r33p-r22p)) + 1i*O13*n21m; 
    n23m_solver.make_step(n23m_t,dt); 
    
    n21p_t = dE21*n21p -1i/2*(pField.*r310 + mField.*r31m) + 1i*O13*n23p;
    n21p_solver.make_step(n21p_t,dt);
    
    n21m_t = dE21*n21m - 1i/2*(mField.*r310 +pField.*r31p) + 1i*O13*n23m;
    n21m_solver.make_step(n21m_t,dt);
    
    if(shb > 0 )
      
        r11p_t = 1i*O13.*(conj(r31m)-r31p) + r33p*W31 + r22p*W21 - G1*r11p;
        r11p_solver.make_step(r11p_t,dt);
        
        r33p_t = 1i*O13.*(r31p-conj(r31m))+1i/2*(mField.*conj(n23p) -conj(pField).*n23m) +  r11p*W13 + r22p*W23 - G3*r33p;
        r33p_solver.make_step(r33p_t,dt);
        
        r22p_t = -1i/2*(mField.*conj(n23p) -conj(pField).*n23m) +  r11p*W12 + r33p*W32 - G2*r22p;
        r22p_solver.make_step(r22p_t,dt);    
        
        r31p_t = dE31*r31p - 1i*O13*(r11p-r33p) -1i/2*(conj(pField).*n21m);
        r31p_solver.make_step(r31p_t,dt);
        r31m_t = dE31*r31m  -1i*O13*conj(r11p-r33p) - 1i/2*(conj(mField).*n21p);
        r31m_solver.make_step(r31m_t,dt);
    end
   
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%
  
    k = -c*l.*ones(N,1);
    pField = forward_solver.make_step(1i*c*n23p,1i*c*n23p_t,k,dt);    
    mField = backward_solver.make_step(1i*c*n23m,1i*c*n23m_t,k,dt);
    
    %set the boundaries... and obtain the final solution
    pField = forward_solver.set_bdry(mField(1),'no');
    mField = backward_solver.set_bdry('no',pField(end));
    
    
    
    r110 = r110_solver.get_latest_solution(); 
    r330 = r330_solver.get_latest_solution(); 
    r220 = r220_solver.get_latest_solution(); 
    
    r310 = r310_solver.get_latest_solution(); 
    n23p = n23p_solver.get_latest_solution(); n23m = n23m_solver.get_latest_solution();
    n21p = n21p_solver.get_latest_solution(); n21m = n21m_solver.get_latest_solution();
    
    if(shb > 0)
        
        r11p = r11p_solver.get_latest_solution();
        r33p = r33p_solver.get_latest_solution();
        r22p = r22p_solver.get_latest_solution();
        r31p = r31p_solver.get_latest_solution(); 
        r31m = r31m_solver.get_latest_solution();
        
    end
    
    t = t+dt;
    iter_ctr = iter_ctr + 1;

end
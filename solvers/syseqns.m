function [ M_prime ] = syseqns( t,M )
    
%     U = M(1,:);
%     V = M(2,:);
%     N_P = M(3,:);
%     N_M = M(4,:); 
%     D_0 = M(5,:); 
%     D_2 = M(6,:);

    global l_0; global alpha;global c; global D_N; global T_1 ; global T_2; global idx_short; global diffusion;
    global shb; global sat_abs; global lambda_L; global p, global m, global f_R,global lTH; 
    n = round(length(M)/6);
    M = reshape(M,[n 6]);   
    M_prime = M;
    
    U = M(:,1);
    V = M(:,2);
    N_P = M(:,3);
    N_M = M(:,4); 
    D_0 = M(:,5); 
    D_2 = M(:,6);
    l = l_0;
    if(sat_abs > 0) 
        l = (l_0 - alpha*(M(:,1)*M(:,1)' + M(:,2)*M(:,2)'));
    end

    M_prime(:,1) = -c*D_N*M(:,1) -1i*c*M(:,3) -c*l*M(:,1);
    M_prime(:,2) = +c*D_N*M(:,2) -1i*c*M(:,4) -c*l*M(:,2);
    M_prime(:,3) = 0.5*1i*(M(:,5).*M(:,1) + conj(M(:,6)).*M(:,2))-...
        M(:,3)/T_2;
    M_prime(:,4) = 0.5*1i*(M(:,5).*M(:,2) + M(:,6).*M(:,1))-...
        M(:,4)/T_2;
    
    g_current =  conj(M(:,1)).*M(:,3) +  conj(M(:,2)).*M(:,4);
    M_prime(1:idx_short,5) = lambda_S(p,m,f_R,t,lTH) +...
        1i*(g_current(1:idx_short) - conj(g_current(1:idx_short))) - M(1:idx_short,5)/T_1;
    M_prime(idx_short+1:end,5) = lambda_L +...
       1i*(g_current(idx_short+1:end) - conj(g_current(idx_short+1:end))) - M(idx_short+1:end,5)/T_1;
   
   if(shb >0)
        M_prime(:,6) = 1i*(conj(M(:,1)).*M(:,4) - conj(M(:,3)).*M(:,2)) - (1/T_1 + diffusion)*M(:,6);
   end
    M_prime = reshape(M_prime,[n*6 1]);

end


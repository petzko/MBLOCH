function [ F ] = constituent_relations_FSA_SG( x, params)
    % x(1) = alpha 
    % x(2) = I0
    % x(3) = tau
    
    T1g = params.T1_g; %ps 
    T2g = params.T2_g; %ps
    WG_sat = params.WG_sat; %ps^-2;
    g0 = params.g0; 
    gamma = params.gamma;
    p = params.p; 
    d_th = params.d_th;

    alpha = 0; %x(1);
    I0 = x(2);
    tau = x(3);
    
    
    W1 = I0*tau/(T1g*WG_sat);
    
    F(1) = alpha-g0*W1+g0*W1^2;
    F(2) = -1/2*g0*W1^2-2*T2g^2/(tau^2)+gamma*I0;
    F(3) = (p-1)*d_th + g0*(-W1+W1^2) +T2g^2/(tau^2);
    

end


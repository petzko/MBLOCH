function [F] = pulse_area_steady_state(X,s,kappa)

% F(1) =  G2 - (G1-g0).*exp(-T/T1g)-g0;
% F(2) =  Q2 - (Q1-q0).*exp(-T/T1a)-q0;
% F(3) =  G1 - G2 - log((exp(G1)-1)./(exp(G2)-1))-DP;
% F(4) =  Q1 - Q2 - log((exp(Q1)-1)./(exp(Q2)-1))-s*log((1+exp(DP))/2);   
% F(5)= s*DP-kappa*(log(1+(1+exp(DP))^s)- log(1+2^s)); 

DP = X(1);
F = s*DP-kappa*(log(1+(1+exp(DP))^s)- log(1+2^s)); 

end


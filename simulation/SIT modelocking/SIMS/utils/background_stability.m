function [F] = background_stability(X,s,T,T1g,T1a,g0,q0,kappa)

% F(1) =  G2 - (G1-g0).*exp(-T/T1g)-g0;
% F(2) =  Q2 - (Q1-q0).*exp(-T/T1a)-q0;
% F(3) =  G1 - G2 - log((exp(G1)-1)./(exp(G2)-1))-DP;
% F(4) =  Q1 - Q2 - log((exp(Q1)-1)./(exp(Q2)-1))-s*log((1+exp(DP))/2);   
% F(5)= s*DP-kappa*(log(1+(1+exp(DP))^s)- log(1+2^s)); 

G1 = X(1);
G2 = X(2);
Q1 = X(3);
Q2 = X(4);
DP = X(5); 

% penalize via the fermi function NEGATIVE G1 and G2 and POSITIVE Q1 and Q2


% Temp = mu/100;
% fermi_G1 = 1/(exp(G1/T)+1); 



% F(1) =  G2 - (G1-g0).*exp(-T/T1g)-g0+0*((abs(abs(G1)-abs(G2)))<tol);
% F(2) =  Q2 - (Q1-q0).*exp(-T/T1a)-q0+0*((abs(abs(Q1)-abs(Q2)))<tol);
F(1) =  G2 - (G1-g0).*exp(-T/T1g)-g0;
F(2) =  Q2 - (Q1-q0).*exp(-T/T1a)-q0;

F(3) =  G1 - G2 - log((exp(G1)-1)./(exp(G2)-1))-DP;
F(4) =  Q1 - Q2 - log((exp(Q1)-1)./(exp(Q2)-1))-s*log((1+exp(DP))/2);   
F(5)= s*DP-kappa*(log(1+(1+exp(DP))^s)- log(1+2^s)); 


% penalize via the fermi function very small DP and negative DP
mu = 1e-4; Temp = mu/2; 
a = sign((1/s-1/kappa)*DP +(s-1)/s*log((1+exp(DP))/2)+(Q1-Q2)/s+(G1-G2)/s);
fermi_DP_a = 1/(exp((DP-mu)/Temp)+1);

% b = sign(s*DP-kappa*(log(1+(1+exp(DP))^s)- log(1+2^s))); 
% fermi_DP_b = b*1/(exp((DP-mu)/Temp)+1);

% F(5)= s*DP-kappa*log((exp(Q2)-1)/(exp(Q1)-1)) + fermi_DP_b; 
F(5)= (1/s-1/kappa)*DP +(s-1)/s*log((1+exp(DP))/2)+(Q1-Q2)/s+(G1-G2)/s +0*fermi_DP_a;
end


function [F] = background_stability_II(X,s,T,g0,q0,Gamma,kappa)

% equations as given in Phys. Rev. A 72, 033808 (2005)
% F(1)= G1-G2*e^(-Gamma * T) - g0/Gamma*(1-e^(-Gamma*T))
% F(2)= Q1-Q2*e^(-T) - q0*(1-e^(-T))
% F(3)= G2 + log(1-(1-e^(-G1))/(e^(-Q1)*(e^(s*DP)-1)+1)^(1/s))
% F(4)= Q2 - log(1+e^(-s*DP)*(e^Q1-1));
% F(5)= DP-kappa*ln((e^(G1)-1)/(e^(G2)-1))

E=exp(1); 

G1 = X(1);
G2 = X(2);
Q1 = X(3);
Q2 = X(4);
DP = X(5); 


F(1)= G1-G2*E^(-Gamma * T) - g0/Gamma*(1-E^(-Gamma*T));
F(2)= Q1-Q2*E^(-T) - q0*(1-E^(-T));
F(3)= G2 + log(1-(1-E^(-G1))/(E^(-Q1)*(E^(s*DP)-1)+1)^(1/s));
F(4)= Q2 - log(1+E^(-s*DP)*(E^Q1-1));
F(5)= DP-kappa*log((E^(G1)-1)/(E^(G2)-1));
% F(5)= DP-kappa*(log(1-(1+E^(s*DP))^1/s)-log(1-2^1/s));
end


function [F] = background_stability_III(X,s,T,g0,q0,kappa,gamma)

G1 = X(1);
G2 = X(2);
Q1 = X(3);
Q2 = X(4);
DP = X(5); 

F(1) = Q2-Q1*exp(-T)-q0*(1-exp(-T));
F(2) = G2-G1*exp(-T*gamma)-g0/gamma*(1-exp(-T*gamma));

F(3) = Q1-Q2-log((exp(Q1)-1)/(exp(Q2)-1))-s*DP;
F(4) = G1-G2-log((exp(G1)-1)/(exp(G2)-1))-gamma/s*log( (1+exp(s*DP))/2 );

F(5) = DP-kappa/gamma*log((exp(G2)-1)/(exp(G1)-1));




end


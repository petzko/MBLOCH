clear; clc; 


g0 = 0:0.1:2.0;
q0 = [0:.1:5];

T=1.87;
s = 25; 
Gamma = 1.33*1e-2;
kappa = .1;

res_1 = zeros(length(g0),length(q0)); 
res_2 = zeros(length(g0),length(q0)); 

options = optimoptions('fsolve','Display','none');

X0 = 1/10*[10,10,10,10,10];
tol = 1e-5; 
for m = 1:length(g0)
    m
    for n=1:length(q0)
         
        
      [X,val] = fsolve('background_stability_II',X0,options,s,T,g0(m),q0(n),Gamma,kappa);
      while (sum(abs(imag(X))) > tol) 
        X0 = [rand(1,2), rand(1,2), rand(1)];
       [X,val] = fsolve('background_stability_II',X0,options,s,T,g0(m),q0(n),Gamma,kappa);
       display('Imag sol')
      end
        X = real(X);
        
        G1 = X(1); G2 = X(2); Q1 = X(3); Q2 = X(4); DP =X(5);  
        res_1(m,n) = G1-Q1+log(kappa) < 0;
        res_2(m,n) = G2-Q2+log(kappa) < 0;
    end
end
res = res_1.*res_2; 
[Gs,Qs] = meshgrid(g0,q0); 
surf(Gs,-Qs,res'); 



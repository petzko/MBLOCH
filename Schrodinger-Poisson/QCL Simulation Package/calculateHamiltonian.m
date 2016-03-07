function [Hext,D,S] = calculateHamiltonian(dV,x,Psi,E)

    nlevels = length(E); 
    Hext = zeros(nlevels);
    S = zeros(nlevels); 
    D = zeros(nlevels); 
    
    for i = 1:nlevels
      Psi(:,i) = Psi(:,i)/sqrt(trapz(x,abs(Psi(:,i)).^2));
    end
    
    for i = 1:nlevels; 
        for j = 1:nlevels;
            S(i,j) =  trapz(x,conj(Psi(:,i)).*dV.*Psi(:,j)) ; 
            D(i,j) = trapz(x,conj(Psi(:,i)).*Psi(:,j))*(i==j); 
            Hext(i,j) = E(j)*D(i,j)*(i==j) + S(i,j);           
        end    
    end


end

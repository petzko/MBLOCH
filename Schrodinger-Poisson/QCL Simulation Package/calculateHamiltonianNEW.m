function [H,D,S] = calculateHamiltonianNEW(Psi,E,dV,x)
    
    %this function calculates the hamiltonian of the full system (witout
    %the TB approximation) in the basis of tight-binding functions Psi. For
    %the calculation one needs to obtain the overlap integral D, the shift
    %integral (and also the coupling strength integral) S and also the
    %resulting H matrix, where H is given by the expression H = E*D+S with
    %E the correctly ordered energies. H and D form an eigenvalue equation 
    % H*a = eps*D*a , where eps are the eigenenergies  of the system,
    % and a are the eigen vectors. Thus to obtain the eigenstates of the
    % full hamiltonian we need to find the eigenvalues of :
    %    D^-1Ha = eps*a
    %As inputs the function takes the eigenstates in the TB basis for a
    %structure consisting of 3 repetitions. Those eigenstates span the
    %whole region however, they are the solution of the tight-binding
    %approximation and are thus confined to a single region. The input
    %variables must follow the folloging convention:
    % Psi - an "length(x)" by "number of levels per period" by "3" array, contianing the 
    % TB wavefunctions divided in periods, where the last index specifies the period where the WFs are localized. 
    % E - a "number of levels" by "3" array containing the eigenenergies in
    % each subperiod
    % dV -  a "length(x)" by "3" array giving the difference between the H_TB hamiltonian in 
    % a specific period and the extended full hamiltonian over the whole
    % domain. 
    nlevels = length(E(:,1));
    H = zeros(3*nlevels);  S = zeros(3*nlevels); D = zeros(3*nlevels); 
    
    % 1 <=> L 
    % 2 <=> 0 
    %   and 
    % 3 <=> R
    
    %normalize the WFs;
    for i = 1:nlevels
      Psi(:,i,1) = Psi(:,i,1)/sqrt(trapz(x,abs(Psi(:,i,1)).^2));
      Psi(:,i,2) = Psi(:,i,2)/sqrt(trapz(x,abs(Psi(:,i,2)).^2));
      Psi(:,i,3) = Psi(:,i,3)/sqrt(trapz(x,abs(Psi(:,i,3)).^2));
    end
    
    % Build first row of H,D and S. 
    
    for A = 1:3
        for B = 1:3
            for mu = 1:nlevels
                for nu = 1:nlevels
                D((A-1)*nlevels + mu, (B-1)*nlevels+nu) = trapz(x,conj(Psi(:,mu,A)).*Psi(:,nu,B));
                S((A-1)*nlevels + mu, (B-1)*nlevels+nu) =  trapz(x,conj(Psi(:,mu,A)).*dV(:,B).*Psi(:,nu,B));
                H((A-1)*nlevels + mu, (B-1)*nlevels+nu) = E(nu,B)*D((A-1)*nlevels + mu, (B-1)*nlevels+nu)+S((A-1)*nlevels + mu, (B-1)*nlevels+nu);
                end
            end
        end
    end
    
    
    


end

function touif = CalcLifetime(Psi,E,i,f,Egx,Vx,mx,dx,hbarwlo,T,epsr)

% This function calculates inverse of the phonon scattering rate 
% or the lifetime between two specified energy levels.
%
% USAGE:
%
% touif = CalcLifetime(Psi,E,i,f,Egx,Vx,mx,dx,hbarwlo,status)
%
% INPUT:
%
% Psi - the wavefunction matrix
% E - eigen energy levels (vector)
% i - initial energy level for the transition
% f - final energy level for the transition
% Egx - bandgap energy mesh
% Vx - potential energy mesh
% mx - electron effective mesh
% dx - grid spacing
% hbarwlo - phonon energy
% status - phonon emission or absorption
%
% OUTPUT:
%
% touif - lifetime between energy levels i and f
%
% Author - Muhammad Anisuzzaman Talukder (anisuzzaman@umbc.edu)
%--------------------------------------------------------------
%
% Theory
% Optical phonon scattering rate is given by -
% 1/touif = (m*e2wlo/(2hbar^2.eps.qif)I(dx (I ...
% (dx'Psii(x).Psif(x).exp(-qif(x-x')).Psii(x')Psif(x') [Ref. 1]
% where qif = sqrt(2m*(Eif-hbarwlo)/hbar^2) is the momentum
% exchanged in the transition.
% For phonon absorption
% qif = sqrt(2m*(Eif+hbarwlo)/hbar^2)
% here Eif is the energy difference between initial and final
% levels, i.e., Eif = Ei-Ef
% For Eif below hbarwlo and at low temperatures, optical phonon
% scattering is forbidden and the lifetime is much longer,
% of the order of hundreds of picoseconds [Ref. 2].

Ei = E(i);
Ef = E(f);
Psii = Psi(:,i);
Psif = Psi(:,f);

nx = length(Psii);

const = Constants();
epsp = const.eps0/(1/epsr.high - 1/epsr.static);
meffi = 0;
mefff = 0;

if Ei > Ef
    status = 'Emission';
else
    status = 'Absorption';
end

% To calculate qif, m* is needed. Here m* is the average or
% expected value in the level. So <m*> is determined.
% <m*> = integrate(Psii.m*.Psii)

if Ei-Ef > hbarwlo

    No = 1/(exp(hbarwlo*const.eV/(const.k*T)) - 1);

    switch status
        case 'Emission'
            
            %for kk = 1 : nx
            %    meffi = meffi + mx(kk) * (Ei-hbarwlo-Vx(kk))/Egx(kk)*dx*Psii(kk)^2;
                %mefff = mefff + mx(kk) * (Ef-hbarwlo-Vx(kk))/Egx(kk)*dx*Psif(kk)^2;
            %end  
            
            meffi = min(mx);
            Eg = min(Egx);
            meffi = meffi*Ei/Eg;
            
            qif = sqrt(2*meffi*const.eV*(Ei-Ef-hbarwlo)/const.hbar^2);
            
            I = 0;
            for kk = 1 : nx
                for jj = 1 : nx
                    I = I + Psii(jj)*Psif(jj)*exp(-qif*abs(kk-jj)*dx)*...
                        Psii(kk)*Psif(kk)*dx*dx;
                end
            end

            invtouif = (No+1)*meffi*abs(const.qe)^3*hbarwlo*I/(2*const.hbar^3*epsp*qif);
      
        case 'Absorption'
            %for kk = 1 : nx
            %    meffi = meffi + mx(kk) * (Ei+hbarwlo-Vx(kk))/Egx(kk)*dx*Psii(kk)^2;
                %mefff = mefff + mx(kk) * (Ef+hbarwlo-Vx(kk))/Egx(kk)*dx*Psif(kk)^2;
            %end
            
            meffi = min(mx);
            Eg = min(Egx);
            meffi = meffi*Ei/Eg;
            
            qif = sqrt(2*meffi*const.eV*(Ei-Ef+hbarwlo)/const.hbar^2);
            
            I = 0;
            for kk = 1 : nx
                for jj = 1 : nx
                    I = I + Psii(jj)*Psif(jj)*exp(-qif*abs(kk-jj)*dx)*...
                        Psii(kk)*Psif(kk)*dx*dx;
                end
            end
        
            invtouif = No*meffi*abs(const.qe)^3*hbarwlo*I/(2*const.hbar^3*epsp*qif);
    end % end switch status
%invtouif
% touif = 1e12/invtouif;
touif = 1./invtouif
Ei - Ef;
else
    touif = 0*1e8;
    %Ei - Ef
end % end if Ei-Ef > hbarwlo

% Reference
% [1] J. Faist, F. Capasso, C. Sirtori, D. Sivco, and Y. Cho, "Intersubband
% Transitions in Quantum Wells: Physics and Device Applications II," 
% Academic Press, pp. 7-8, 2000.
% [2] J. Faist, C. Sirtori, F. Capasso, L. Pfeiffer, and K. W. West,
% Appl. Phys. Lett. 64, 872, 1994.
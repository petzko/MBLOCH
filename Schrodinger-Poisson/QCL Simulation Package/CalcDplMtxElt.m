function dplmtxelt = CalcDplMtxElt(Ei,Ef,Psii,Psif,Egx,Vx,mx)

% This function calculates the dipole matrix element
% between two specified energy levels
%
% USAGE:
%
% dplmtxelt = CalcDplMtxElt(Ei,Ef,Psii,Psif,Egx,Vx,mx)
%
% INPUT:
%
% Ei - initial energy level
% Ef - final energy level
% Psii - wavefunction of the energy level i
% Psif - wavefunction of the energy level f
% Egx - bandgap energy mesh
% Vx - potential energy mesh
% mx - electron effective mass mesh
%
% OUTPUT:
%
% dplmtxelt - dipole matrix element beween energy levels i and f
%
% Author - Muhammad Anisuzzaman Talukder (anisuzzaman@umbc.edu)
%--------------------------------------------------------------
%
% Theory
%
% Generally, the definition of the dipole matrix element between
% state i and f is
% zif = < Psii|z|Psif >
% However, since we solve a one-dimensional Schrodinger equation,
% which includes the energy dependent effective mass
% m*(E) = m*(E=0)(1+(E-V)/Eg),
% the wavefunctions that we compute with this approach are
% not orthogonal, which makes the definition of the dipole matrix
% element between the states i and f somewhat problematic [Ref. 1].
%
% The solution to this problem is to go back to the two-band
% model and compute the matrix element including the valence
% band part. The dipole matrix element now reads [Ref. 2]
% zif = hbar/(2(Ef-Ei) < Psii|pz(1/m*(Ei,z)) + (1/m*(Ef,z))pz|Psif >
% where the momentum operator pz is defined, as usual, 
% pz = -ihbar(del/delz).
% To account for the underlying valence band component, 
% the wave functions Psii and Psif must be normalized according to
% 1 = < Psii|1 + (E-V(z))/(E-V(z)+Eg(z))|Psif >

const = Constants();
nx = length(Psii);
wif = (Ei-Ef)*const.eV/const.hbar;

% Following section of code (for loop) includes the 
% nonparabolicity due to the two energy levels Ei and Ef
% in the mass mesh. [Ref. 1,3]

for ii = 1 : nx
    mxi(ii) = mx(ii) * (1 + (Ei-Vx(ii))/Egx(ii));
    mxf(ii) = mx(ii) * (1 + (Ef-Vx(ii))/Egx(ii));
end

for kk = 1 : nx-1
    dpl(kk) = 0.5*(const.hbar/(mxi(kk)*wif)*(Psii(kk) + Psii(kk+1))/2*(Psif(kk+1) - Psif(kk))+...
                const.hbar/(mxf(kk)*wif)*(Psii(kk) + Psii(kk+1))/2*(Psif(kk+1) - Psif(kk)));
end

dplmtxelt = sum(dpl);

% Reference
% [1] J. Faist, F. Capasso, C. Sirtori, D. Sivco, and A. Cho,
% "Intersubband Transitions in Quantum Wells: Physics and Device
% Applications II," Academic Press, Ch. 1, p. 6-7, 2000.
% [2] C. Sirtori, F. Capasso, and J. Faist, "Nonparabolicity and
% a sum rule associated with bound-to-bound and bound-to-continuum
% intersubband transitions in quantum wells," Phys. rev. B, 50, 
% p. 8663, 1994.
% [3] D. Nelson, R. Miller, and D. Kleinman, "Band nonparabolicity
% effects in semiconductor quantum wells," Phsy. Rev. B 35, 7770, 1987.
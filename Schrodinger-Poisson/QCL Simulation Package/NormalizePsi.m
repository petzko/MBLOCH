function Psi = NormalizePsi(Psi,E,Egx,Vx,dx,type)

% This function normalizes the wavefunction
%
% USAGE:
%
% NormalizePsi(Psi,E,nlevel,Egx,Vx,dx)
%
% INPUT:
%
% Psi - wavefunction matrix
% E - eigen energy levels (vector)
% Egx - bandgap energy mesh
% Vx - potential energy mesh
% mx - electron effective mass mesh
%
% OUTPUT:
%
% Psi - wavefunction matrix

% Author - Muhammad Anisuzzaman Talukder (anisuzzaman@umbc.edu)
%--------------------------------------------------------------

nlevel = length(E);
nx = length(Psi(:,1));

switch type
    case 1
        for kk = 1 : nlevel
            Psikk = Psi(:,kk);
            normfac = 0;
            for ii = 1 : nx
                normfac = normfac + Psikk(ii)*Psikk(ii);
            end
            Psi(:,kk) = (1/sqrt(dx*normfac))*Psikk;
        end
    case 2
        for kk = 1 : nlevel
            Psikk = Psi(:,kk);
            normfac = 0;
            for ii = 1 : nx
                normfac = normfac + Psikk(ii)*(1+((E(kk)-Vx(ii))/(E(kk)-Vx(ii)+Egx(ii))))*Psikk(ii);
            end
    
            Psi(:,kk) = (1/sqrt(dx*normfac))*Psikk;
        end
end % end switch type
    
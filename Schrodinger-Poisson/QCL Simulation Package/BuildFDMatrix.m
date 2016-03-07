function A = BuildFDMatrix(dx,Egx,Vx,Vxflat,mx,boundary,E)

% This function builds the finite difference matrix for the
% quantum well. Using the potential mesh specified in Vx
% and the boundary conditions specified in 'boundary', this
% function constructs a sparse matrix A representing the 
% differential operator for the eigenvalues.
%
% USAGE:
%
% A = BuildMatrix(dx,Egx,Vx,mx,boundary,E)
%
% INPUT:
%
% dx - grid size
% Egx - bandgap potential mesh
% Vx - potential mesh
% mx - mass mesh
% boundary - boundary condition to be applied at the edges of 
%            the computaion window.
% The following boundary conditions are supported:
%   'zero' - wavefunction is zero immediately outside of the 
%            boundary
%   'periodic' - wavefunction is periodic outside of the boundary
% E - energy of the state
%
% OUTPUT:
%
% A - sparse matrix representing differential operator for the
%     eigenvalue problem
%
% Author - Muhammad Anisuzzaman Talukder (anisuzzaman@umbc.edu)

const = Constants();
nx = length(Vx);

if nargin == 7
    for kk = 1 : nx
        mx(kk) = mx(kk) * (1 + (E-Vx(kk))/Egx(kk));
    end
end

aw = ones(nx-1,1)/dx^2;
ae = ones(nx-1,1)/dx^2;
ap = -ones(nx,1)/dx^2;

mxn = zeros(1,nx);

for kk = 1 : nx
    if kk == nx
        mxn(kk) = (mx(kk) + mx(1))/2;
    else
        mxn(kk) = (mx(kk) + mx(kk+1))/2;
    end
end

for ii = 1 : nx-1
    aw(ii) = aw(ii)/mxn(ii);
    ae(ii) = ae(ii)/mxn(ii);
end

aw = [aw ; 0];
ae = [0 ; ae];

for ii = 1 : nx
    if ii == 1
        ap(ii) = ap(ii)*(1/mxn(ii)+1/mxn(nx));
    else
        ap(ii) = ap(ii)*(1/mxn(ii)+1/mxn(ii-1));
    end

end
ap = (-const.eV*2/const.hbar^2)*Vx' + ap;

switch boundary
    case 'zero'
        A = spdiags([aw,ap,ae],[-1 0 1],nx,nx);
    case 'periodic'
        al = 1/(mxn(1)*dx^2);
        ar = 1/(mxn(nx)*dx^2);
        A = spdiags([[al;zeros(nx-1,1)],aw,ap,ae,[zeros(nx-1,1);ar]],[-(nx-1) -1 0 1 (nx-1)],nx,nx);
    otherwise
        error('Boundary condition is not set properly');
end
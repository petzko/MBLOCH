function [Psi,E] = SolveEigWaves(A,nx,nlevel,guess)

% This function calculates the eigenvalues/energy states and
% the wavefunctions. It uses MATLAB function eigs to calculate
% few eigenvalues/energy states of the matrix A.
%
% USAGE:
%
% [Psi,E] = SolveEigWaves(A,nx)
%
% INPUT:
%
% A - sparse matrix containing the finite-difference
%     representation of the differential operator generated by
%     BuildMatrix.m
% nx - dimension of the finite differnce mesh
% nlevels - number of eigenvalues/energy states
% guess - scalar shift to apply when calculating the eigenvalues. 
%         This routine will return the eigenpairs which are
%         closest to this guess in magnitude
%
% OUTPUT:
%
% Psi - 1 D electron wavefunction for each calculated eigenvalue
% E - vector of the eigenvalues/energy states
%
% Author - Muhammad Anisuzzaman Talukder (anisuzzaman@umbc.edu)

const = Constants();

options.tol = 1E-12;
options.disp = 0;           % suppress output

%note that the WFs are accurate up to a random phase!
[v,d] = eigs(A,speye(size(A)),nlevel,guess,options);

E = diag(d)/(2/const.hbar^2);  % energy levels in joules
E = -E/const.eV;           % energy levels in electron-volts or volts

Psi = zeros(nx,nlevel);

Psi = v;
function [structure,dipoles,H] = sbs(outdir,efield,infile,nlevel,T,nonpar_flag)
% sbs solves the energy states and the wave-functions
% of quantum cascade laser
% sbs - solve band structure
%--------------------------------------------------------------------------
%
% USAGE:
%
% function sbs(efield,infile,nlevel,T,nonpar_flag)
%
% INPUT:
%
% efield - applied external electric field
%   efield.value - value of applied external electric field
%   efield.direction - direction of applied external electric field
% infile - input file for the structure of the QCL
% nlevel - number of states to be solved
% T - temperature
% nonpar_flag - non-parabolicity flag
%
% Author: Muhammad Anisuzzaman Talukder (anisuzzaman@umbc.edu) 
% -------------------------------------------------------------------------

% Simulation parameters
%----------------------
dx = 1e-10;                 % grid size in m


%efield.direction = '<-';    % polarity of the external applied electric potential
boundary = 'zero';          % boundary condition
%  boundary = 'periodic';
%nlevel = 20;                % number of eigenvalues/energy states
guess = 1e-6;               % eigenvalues will be calculated closest to this value

ilife = 3;                  % initial energy level for lifetime calculation
flife = 1;                  % final energy level for lifetime calculation

%T = 300% temperature, K

% purge parameters
purge.flag = 0;             % determines whether purging effect will be included
purge.time = [0 0];         % purge.time(1) - after InAlAs; purge.time(2) - after InGaAs
purge.npdepth = 30e-10;     % nonpurged penetration at the interface
purge.maxtime = 20;         % maximum purging time required for the sharp interface


errorcontrol_flag = 0;

% Following flags control the calculation of
% different parameters for the QCL.
% nonpar_flag = 1;
dipole_flag = 1;
neff_flag = 0;
lifetime_flag = 0;
dope_flag = 0;

% Following flags control plotQCL function
wavefunc_flag = 0;
modulisqr_flag = 1;
modulisqrtrunc_flag = 0;

%format 'long';

fprintf(1,'Simulation Starts\n');

fid = fopen(strcat(outdir,'/','SolveBandStructure.txt'),'w');
fprintf(fid,'%s\n','Modeling Quantum Cascade Laser');
fprintf(fid,'%s\n','------------------------------');
fprintf(fid,'\n%s %s\n','Input file:',infile);
fprintf(fid,'\n%s\n','Operating parameters');
fprintf(fid,'%s\n','--------------------');
fprintf(fid,'%s\n','Applied electric field');
fprintf(fid,'    %s %2.2f %s\n','value:',efield.value,'kV/cm');
if efield.value ~= 0
    fprintf(fid,'    %s %s\n','directions:',efield.direction);
end
fprintf(fid,'%s %3.0f %s\n','Temperature:',T,'K');
fclose(fid);

fid = fopen(infile, 'r');
structure = fscanf(fid, '%g ', [4 inf]);    
fclose(fid);

[x,nx,Egx,Vx,mx,dopx,rindexx,epsr,hbarwlo,Vxflat] = QCLMesh(outdir,structure,dx,efield,T,purge);
A = BuildFDMatrix(dx,Egx,Vx,Vxflat,mx,boundary);
[Psi,E] = SolveEigWaves(A,nx,nlevel,guess);

if errorcontrol_flag == 1
    error = 1;
    while error > 1e-2
        dx = dx/2;
        Em = E;
        [x,nx,Egx,Vx,mx,dopx,rindexx,epsr,hbarwlo,Vxflat] = QCLMesh(outdir,structure,dx,efield,T,purge);
        A = BuildFDMatrix(dx,Egx,Vx,Vxflat,mx,boundary);
        [Psi,E] = SolveEigWaves(A,nx,nlevel,guess);
        error = abs(max(E-Em));
    end
end

nx

if nonpar_flag == 1

for kk = 1 : nlevel
    fixing=kk;
    En = E(kk);
    err = 1;
    eps = 1e-4;
   
    while err > eps 

        A = BuildFDMatrix(dx,Egx,Vx,Vxflat,mx,boundary,En);
        [Psitemp,Etemp] = SolveEigWaves(A,nx,nlevel,guess);

        err = abs(En - Etemp(kk));
        
        if err < eps
            Psi(:,kk) = Psitemp(:,kk);
        end

        En = Etemp(kk);
                        
    end % end while err > eps
    
    E(kk) = En;
    
end
end % end if nonpar == 1

% If dope_flag is set to 1, self-consistent solution
% for the quantum structure is determined for the 
% modulation doping solving Schrodinger equation and 
% Poission's equation iteratively until the values conerge.
Vx_dope = zeros(1,nx);
Vx_dopeold = zeros(1,nx);
Vxold = Vx;

if dope_flag == 1
    Psi = NormalizePsi(Psi,E,Egx,Vx,dx,1);
    Vx_dope = CalcVx_dope(Psi,E,dopx,rindexx,dx);
    Vx = Vxold + Vx_dope;
    convergence = abs(sum(Vx_dope)-sum(Vx_dopeold));

    for abc = 1 : 1
    %while convergence > 1e-3
    
        A = BuildFDMatrix(dx,Egx,Vx,Vxflat,mx,boundary);
        [Psi,E] = SolveEigWaves(A,nx,nlevel,guess);
     
        if nonpar == 1
            for kk = 1 : nlevel
                En = E(kk);
                err = 1;
                eps = 1e-2;

                while err > eps

                    A = BuildFDMatrix(dx,Egx,Vx,Vxflat,mx,boundary,En);
                    [Psitemp,Etemp] = SolveEigWaves(A,nx,nlevel,guess);

                    err = abs(En - Etemp(kk));
        
                    if err < 1e-2
                        Psi(:,kk) = Psitemp(:,kk);
                    end

                    En = Etemp(kk);
                        
                end % end while err > eps
    
                E(kk) = En;

            end % end for kk = 1 : nlevel
        end % end if nonpar == 1

        Psi = NormalizePsi(Psi,E,Egx,Vx,dx,1);
        Vx_dopeold = Vx_dope;
        Vx_dope = CalcVx_dope(Psi,E,dopx,rindexx,dx);
        
        Vx = Vxold + Vx_dope;
        
        convergence = abs(sum(Vx_dope)-sum(Vx_dopeold));
        convergence
    
    end % end while convergence > 1e-10
end % end if dope_flag == 1

%E

dipoles = zeros(nlevel);
if dipole_flag == 1
      Psi = NormalizePsi(Psi,E,Egx,Vx,dx,2);
    for idpl = 1:nlevel
        for fdpl = idpl:nlevel 
            dipoles(idpl,fdpl) = 1e10*abs(CalcDplMtxElt(E(idpl),E(fdpl),Psi(:,idpl),Psi(:,fdpl),Egx,Vx,mx));
            dipoles(fdpl,idpl) = dipoles(idpl,fdpl);
        end
    end
    
end

% Effective refractive index for the quantum structure
% is calculated if neff_flag is set to 1.
if neff_flag == 1
    neff = Calcneff(rindexx,dx);
end

% Phonon scattering rate is calculated between two 
% energy eigenvalues i and f is lifetime_flag is set to 1.
if lifetime_flag == 1
life_times = zeros(nlevel); 
% 
%    Psi = NormalizePsi(Psi,E,Egx,Vx,dx,1);
%    for i_idx = 1:nlevel
%        for j_idx = 1:nlevel
%         touifemit = CalcLifetime(Psi,E,i_idx,j_idx,Egx,Vx,mx,dx,hbarwlo,T,epsr);
%         life_times(i_idx,j_idx) = real(touifemit);
%         touif = touifemit
%        end
%    end

% life_times
% scat_rates = 1./life_times;
% scat_rates
end
%plotQCL(Psi,E,Vx,x,wavefunc_flag,modulisqr_flag,modulisqrtrunc_flag);

fid = fopen(strcat(outdir,'/','SolveBandStructure.txt'),'a');
fprintf(fid,'\n%s\n','Calculated values');
fprintf(fid,'%s\n','-----------------');
fprintf(fid,'%s\n','Energy Eigenvalues:');
fprintf(fid,'    %2.6f\n',E);
if dipole_flag == 1
    fprintf(fid,'\n%s\n','Dipole matrix element');
    fprintf(fid,'    %s %2.0f %2.0f\n','Levels:',idpl,fdpl);
    fprintf(fid,'    %s %2.4f %s\n','Ei-Ef =',E(idpl)-E(fdpl),'eV');
end
if lifetime_flag == 1
    fprintf(fid,'\n%s\n','Lifetime');
    fprintf(fid,'    %s %2.0f %2.0f\n','Levels:',ilife,flife);
    fprintf(fid,'    %s %2.4f %s\n','Ei-Ef =',E(ilife)-E(flife),'eV');
    fprintf(fid,'    %s %3.6f %s\n','touif =',touif,'ps');
end
if neff_flag == 1
    fprintf(fid,'\n%s\n','Effective refractive index');
    fprintf(fid,'    %s %3.4f\n','neff =',neff); 
end
fclose(fid);

fid = fopen(strcat(outdir,'/','Psiout.txt'),'w');
fprintf(fid,'%2.2e\n',Psi);
fclose(fid);

fid = fopen(strcat(outdir,'/','EnergyValues.txt'),'w');
fprintf(fid,'%4.6f\n',E);
fclose(fid);

fid = fopen(strcat(outdir,'/','Vxmesh.txt'),'w');
fprintf(fid,'%4.6f\n',Vx);
fclose(fid);

fid = fopen(strcat(outdir,'/','Vxflat.txt'),'w');
fprintf(fid,'%4.6f\n',Vxflat);
fclose(fid);

fid = fopen(strcat(outdir,'/','mxmesh.txt'),'w');
fprintf(fid,'%4.6e\n',mx);
fclose(fid);

fid = fopen(strcat(outdir,'/','Egxmesh.txt'),'w');
fprintf(fid,'%4.6f\n',Egx);
fclose(fid);

fid = fopen(strcat(outdir,'/','dopxmesh.txt'),'w');
fprintf(fid,'%4.6e\n',dopx);
fclose(fid);

fid = fopen(strcat(outdir,'/','rindexxmesh.txt'),'w');
fprintf(fid,'%4.6f\n',rindexx);
fclose(fid);

fid = fopen(strcat(outdir,'/','epsr.txt'),'w');
fprintf(fid,'%4.6f\n',epsr.high);
fprintf(fid,'%4.6f\n',epsr.static);
fclose(fid);

fid = fopen(strcat(outdir,'/','hbarwlo.txt'),'w');
fprintf(fid,'%4.6f\n',hbarwlo);
fclose(fid);

fid = fopen(strcat(outdir,'/','Temp.txt'),'w');
fprintf(fid,'%4.6f\n',T);
fclose(fid);

fid = fopen(strcat(outdir,'/','GridSize.txt'),'w');
fprintf(fid,'%4.6e\n',dx);
fclose(fid);

fid = fopen(strcat(outdir,'/','x.txt'),'w');
fprintf(fid,'%4.6e\n',x);
fclose(fid);

fid = fopen(strcat(outdir,'/','Infile.txt'),'w');
fprintf(fid,'%s',infile);
fclose(fid);
H = A;
% This progrma will plot the conduction band edge potential and the solved
% moduli-squared wavefunctions for a quantum cascad lasers (QCLs). The QCL
% strucute is in a file in the folder "Input." Solved wavefunctions, energy
% levels, and other quantitites are saved in the folder "Output."
% This QCL simulation tool has been written by 
% Muhammad A. Talukder, Dept. of CSEE, UMBC
% email: anisuzzaman@umbc.edu
%--------------------------------------------------------------------------

% clear all;
% close all;

efield.value = 10.3;                  % electric field value
efield.direction = '->';            % electric field direction
indir = 'Input';
infile = strcat(indir,'/','QCL183S.m');  % input file
outdir = 'Output';
nlevel = 10;                        % number of states to be solved
T = 50;                            % Temperature
nonpar_flag = 1;                    % non-parabolicity flag

%--------------------------------------------------------------------------
sbs(outdir,efield,infile,nlevel,T,nonpar_flag);

A = load(strcat(outdir,'/','Psiout.txt'));
E = load(strcat(outdir,'/','EnergyValues.txt'));
Vx = load(strcat(outdir,'/','Vxmesh.txt'));
Egx = load(strcat(outdir,'/','Egxmesh.txt'));
mx = load(strcat(outdir,'/','mxmesh.txt'));
ep = load(strcat(outdir,'/','epsr.txt'));
hbarwlo = load(strcat(outdir,'/','hbarwlo.txt'));
T = load(strcat(outdir,'/','Temp.txt'));
dx = load(strcat(outdir,'/','GridSize.txt'));
x = load(strcat(outdir,'/','x.txt'));
nx = length(Vx);
const = Constants();

Psi = zeros(nx,nlevel);

for kk = 1 : nlevel
    sindex = (kk-1)*nx + 1;
    findex = sindex+nx-1;
    Psi(:,kk) = A(sindex : findex);
end
% f_1_name = ['WF_' num2str(efield.value) '.txt'];
% f_2_name = ['POT_' num2str(efield.value) '.txt'];
% plotQCL(Psi,E,Vx,x,0,1,0);


efield.value = efield.value * 1e3 * 1e2;

switch efield.direction
    case '->'
        vx = (1:1:nx)*dx*efield.value;
    case '<-'
        vx = (nx:-1:1)*dx*efield.value;
    otherwise
        error('Direction of the applied electric field is not set properly'); 
end



VB = Vx(1); % in eV! 
epsilon = 1000E-10; % 6 nm
Lp = 590E-10; %% 54nm

Vx_TB = (vx').*(x<Lp)+ (VB+vx').*(x>=Lp).*(x<=Lp+epsilon);
plot(x,Vx,x,Vx_TB);
dV = (Vx_TB-Vx);
plot(x,dV);

% Psi_tb = zeros(nx,nlevel);
% E_tb = zeros(length(E),1); 
% x_Angstrom = x*1E10;
% 
% for jj = 1 : nlevel
%     E_tb(jj) = E(jj) +VB*trapz(dV.*abs(Psi(:,jj)).^2,x_Angstrom);
%     
%     Psi_tb(:,jj) = Psi(:,jj);
%     for ii = 1:nlevel
%         if ii == jj 
%             continue; 
%         else
%             Psi_tb(:,jj) = Psi_tb(:,jj) + VB/(E(jj)-E(ii))*trapz(dV.*conj(Psi(:,ii)).*Psi(:,jj),x_Angstrom)*Psi(:,ii);
%         end
%     end
% end
% % plotQCL(Psi,E,Vx,x,0,1,0,'EB');
% plotQCL(Psi_tb,E,Vx,x,0,1,0,'TB');

% write_WF(x*1E10,Psi,E,Vx,efield.value,f_1_name,f_2_name)
% 
% movefile(f_1_name,[ pwd '\post processing\' f_1_name] )
% movefile(f_2_name,[ pwd '\post processing\' f_2_name] )
% 


% This progrma will plot the conduction band edge potential and the solved
% moduli-squared wavefunctions for a quantum cascad lasers (QCLs). The QCL
% strucute is in a file in the folder "Input." Solved wavefunctions, energy
% levels, and other quantitites are saved in the folder "Output."
% This QCL simulation tool has been written by 
% Muhammad A. Talukder, Dept. of CSEE, UMBC
% email: anisuzzaman@umbc.edu
%--------------------------------------------------------------------------

clear all;
close all;
biases =6.3;
zmn_1 = zeros(length(biases),1);
En_1 = zeros(length(biases),1);
En_2 = zeros(length(biases),1); 
En_3 = zeros(length(biases),1);



ctr = 1; 
for b = biases
    b
    efield.value = b;                  % electric field value
    efield.direction = '->';            % electric field direction
    indir = 'Input';
    infile = strcat(indir,'/','QCL183STB.m');  % input file
    outdir = 'Output';
    nlevel = 3;                        % number of states to be solved
    T = 1;                            % Temperature
    nonpar_flag = 1;                    % non-parabolicity flag
    
    %--------------------------------------------------------------------------
    [struct, dipole ] = sbs(outdir,efield,infile,nlevel,T,nonpar_flag);
    zmn_1(ctr) = dipole(1,2); 
    
    
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
    
    En_1(ctr) = E(1); 
    En_2(ctr) = E(2); 
    En_3(ctr) = E(3); 
    
    ctr = ctr+1;
    
%     f_1_name = ['WF_' num2str(efield.value) '.txt'];
%     f_2_name = ['POT_' num2str(efield.value) '.txt'];
    plotQCL(Psi,E,Vx,x,0,1,0,infile);  

end

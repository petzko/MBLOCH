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
biases = 9:.05:12;
zmn_1 = zeros(length(biases),2);
zmn_2 = zeros(length(biases),2);
zmn_3 = zeros(length(biases),2);


ctr = 1; 
for b = biases
    b
    efield.value = b;                  % electric field value
    efield.direction = '->';            % electric field direction
    indir = 'Input';
    infile = strcat(indir,'/','RT2WELL.m');  % input file
    outdir = 'Output';
    nlevel = 7;                        % number of states to be solved
    T = 50;                            % Temperature
    nonpar_flag = 1;                    % non-parabolicity flag

    %--------------------------------------------------------------------------
    [struct, dipole ] = sbs(outdir,efield,infile,nlevel,T,nonpar_flag);
    zmn_1(ctr,1) = dipole(5,4); 
    zmn_1(ctr,2) = dipole(5,3); 
    zmn_2(ctr,1) = dipole(6,4); 
    zmn_2(ctr,2) = dipole(6,3); 
    zmn_3(ctr,1) = dipole(7,4); 
    zmn_3(ctr,2) = dipole(7,3); 
    ctr = ctr +1;
    
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
    f_1_name = ['WF_' num2str(efield.value) '.txt'];
    f_2_name = ['POT_' num2str(efield.value) '.txt'];
%     plotQCL(Psi,E,Vx,x,0,1,0,infile);  

    write_WF(x*1E10,Psi,E,Vx,efield.value,f_1_name,f_2_name)
    movefile(f_1_name,[ pwd '\post processing\' f_1_name] )
    movefile(f_2_name,[ pwd '\post processing\' f_2_name] ) 
end

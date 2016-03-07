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
bias = 11;
AC_energies_left = zeros(length(bias),2);
AC_energies_right = zeros(length(bias),2);
const = Constants();
nper = 3;
dEnergies = zeros(length(bias),4);

ctr =1;
for b = bias
    b
    efield.direction = '->';            % electric field direction
    efield.value = b;                  % electric field value
    indir = 'Input';
    infile = strcat(indir,'/','QCL2WELL.m');  % input file
    outdir = 'Output';
    nlevel = 7;                        % number of states to be solved
    T = 50;                            % Temperature
    nonpar_flag = 0;                    % non-parabolicity flag
    
    %--------------------------------------------------------------------------
    [structure_ext,dipole_ext,FD_ext ] =  sbs(outdir,efield,infile,nlevel,T,nonpar_flag);
    CB_W = structure_ext(1,1)*1E-10; %coupled basis confinement barrier width
    Lp = (sum(structure_ext(1,1:end-1))*1E-10)/nper;
    
    
    A = load(strcat(outdir,'/','Psiout.txt'));
    E_ext = load(strcat(outdir,'/','EnergyValues.txt'));
    dx = load(strcat(outdir,'/','GridSize.txt'));
    Vx = load(strcat(outdir,'/','Vxmesh.txt'));
    x = load(strcat(outdir,'/','x.txt'));
    nx = length(Vx);
    
    Psi = zeros(nx,nlevel);
    
    for kk = 1 : nlevel
        sindex = (kk-1)*nx + 1;
        findex = sindex+nx-1;
        Psi(:,kk) = A(sindex : findex);
    end
    
    efield.value = efield.value * 1e3 * 1e2;
    
    switch efield.direction
        case '->'
            vx = (1:1:nx)*dx*efield.value;
        case '<-'
            vx = (nx:-1:1)*dx*efield.value;
        otherwise
            error('Direction of the applied electric field is not set properly');
    end
    
    
    % TB Part
    efield.value = b;                  % electric field value
    indir = 'Input';
    
    infile = strcat(indir,'/','QCL2WELL.m');  % input file
    outdir = 'Output';
    nlevel = 3;                        % number of states to be solved
    T = 50;                            % Temperature
    nonpar_flag = 0;                    % non-parabolicity flag
    
    %--------------------------------------------------------------------------
    [ structure_TB,dipole_TB,FD_TB] = sbs(outdir,efield,infile,nlevel,T,nonpar_flag);
    TB_W = structure_TB(1,1)*1E-10;
    %TB Basis confinement barrier width in Angstrom CB-TB barrier width difference;
    dW = CB_W-TB_W; dL =(dW);
    
    
    ATB = load(strcat(outdir,'/','Psiout.txt'));
    E = load(strcat(outdir,'/','EnergyValues.txt'));
    VxTB = load(strcat(outdir,'/','Vxmesh.txt'));
    dxTB = load(strcat(outdir,'/','GridSize.txt'));
    xTB = load(strcat(outdir,'/','x.txt'));
    nxTB = length(VxTB);
    PsiTB_old = zeros(nxTB,nlevel);
    PsiTB = zeros(nx,nper*nlevel);
    ETB = zeros(nper*nlevel,1);
    
    
    
    VpM = b*1E5*Lp; % voltage per module.
    for kk = 1:nlevel
        sindex = (kk-1)*nxTB + 1;
        findex = sindex+nxTB-1;
        baselvl = ATB(sindex : findex);
        PsiTB_old(:,kk)=baselvl;
        baselvl = interp1(xTB+dL,baselvl,x,'linear',0);
        
        PsiTB(:,kk)=baselvl;
        ETB(kk) = E(kk)+b*1E3*1E2*(dL);
        
        %shift the remaining levles
        for i = 1:nper-1
            shiftlvl = interp1(x+i*Lp,baselvl,x,'linear',0);
            PsiTB(:,kk+i*nlevel) = shiftlvl;
            ETB(kk+i*nlevel) = ETB(kk)+i*b*1E3*1E2*Lp;
        end
    end
    
    
    
    %left WFS
    %renormalize wfs;
    for i = 1:nper*nlevel
        Psi_i = PsiTB(:,i);
        PsiTB(:,i) = Psi_i/sqrt(trapz(x,abs(Psi_i).^2));
    end
    
    %x in m and Vx in eV, bias in kV/cm
    VB = Vx(1);
    if efield.direction == '->'
        V0 = Vx-x*b*1E3*1E2;
    else
        V0 = Vx+x*b*1E3*1E2;
    end
    V0TB = V0.*(x>=Lp+CB_W/2).*(x<=2*Lp+CB_W/2);
    dV = V0-V0TB;
    
    %         plotQCL(Psi,E_ext,Vx,x,0,1,0,'ext');
    %         plotQCL(PsiTB,ETB,Vx,x,0,1,0,'TB');
    AC_left = zeros(nlevel); AC_right = zeros(nlevel);
    
    for i =1:nlevel
        Psi_i = PsiTB(:,i);
        for j = 1:nlevel
            Psi_j = PsiTB(:,j+nlevel);
            AC_left(i,j) = trapz(x,conj(Psi_i).*dV.*Psi_j);
        end
    end
    
    for i =1:nlevel
        Psi_i = PsiTB(:,i+nlevel);
        for j = 1:nlevel
            Psi_j = PsiTB(:,j+2*nlevel);
            AC_right(i,j) = trapz(x,conj(Psi_i).*dV.*Psi_j);
        end
    end
    ULL = 3; INJ1= 1; INJ2 = 1; LLL1 = 2; LLL2 = 2; 
    AC_energies_left(ctr,1) = AC_left(ULL,INJ1);  AC_energies_left(ctr,2) = AC_left(ULL,INJ2);
    AC_energies_right(ctr,1) = AC_right(ULL,INJ1);  AC_energies_right(ctr,2) = AC_right(ULL,INJ2);
    dEnergies(ctr,1) = ETB(nlevel+INJ1)-ETB(ULL); 
    dEnergies(ctr,2) = ETB(nlevel+INJ2)-ETB(ULL); 
    dEnergies(ctr,3) = ETB(ULL)-ETB(LLL1); 
    dEnergies(ctr,4) = ETB(ULL)-ETB(LLL2); 
    
    
    ctr = ctr + 1;
    
end


legend_info = {'left ull<-> inj1 ','left ull <-> inj2','right ull<-> inj1 ','right ull <-> inj2'};
plot(bias,abs(AC_energies_left),bias,abs(AC_energies_right)); 
legend(legend_info); 

figure; legend_info = {'dE inj1<-> ull ','dE inj2 <-> ull','dE ull <-> lll1 ','dE ull <-> lll2'};
plot(bias,dEnergies); legend(legend_info); 


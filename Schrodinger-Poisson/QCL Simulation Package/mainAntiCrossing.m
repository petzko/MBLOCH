% This progrma will plot the conduction band edge potential and the solved
% moduli-squared wavefunctions for a quantum cascad lasers (QCLs). The QCL
% strucute is in a file in the folder "Input." Solved wavefunctions, energy
% levels, and other quantitites are saved in the folder "Output."
% This QCL simulation tool has been written by 
% Muhammad A. Talukder, Dept. of CSEE, UMBC
% email: anisuzzaman@umbc.edu
%--------------------------------------------------------------------------

clear all;
% close all;
bias = 9:0.1:12;
% titles = cell(length(bias),1);
nlevel = 7;   % number of states to be solved
% Mov = zeros(length(bias),1); 

Energies = zeros(length(bias),nlevel);
ctr = 1;
for b = bias
b
efield.value = b;                  % electric field value
efield.direction = '->';            % electric field direction
indir = 'Input';
infile = strcat(indir,'/','QCL183S.m');  
outdir = 'Output';
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
    

    inj = [7,6,5]; 
    lll = [4,3];
    A = zeros(length(inj),length(lll));
    for i = 1:length(inj)
        i_1 = inj(i); 
        j_1 = lll(1); j_2 = lll(2); 
        A(i,1) = CalcDplMtxElt(E(i_1),E(j_1),Psi(:,i_1),Psi(:,j_1),Egx,Vx,mx);
        A(i,2) = CalcDplMtxElt(E(i_1),E(j_2),Psi(:,i_1),Psi(:,j_2),Egx,Vx,mx);
    end

    [r,c] = find(abs(A)== max(abs(A(:))));
    ull = inj(r);
    injectors = setdiff(inj,ull);
    l_1 = lll(c);
    l_2 = setdiff(lll,l_1);

    Psi_new_small = zeros(nx,5); E_small = zeros(5,1);
    Psi_new_small(:,1) = Psi(:,injectors(1)); E_small(1) = E(injectors(1));
    Psi_new_small(:,2) = Psi(:,injectors(2)); E_small(2) = E(injectors(2));
    Psi_new_small(:,3) = Psi(:,ull); E_small(3) = E(ull); 
    Psi_new_small(:,4) = Psi(:,l_1); E_small(4) = E(l_1); 
    Psi_new_small(:,5) = Psi(:,l_2); E_small(5) = E(l_2);

    dEnergies(ctr,1) = E_small(1)-E_small(3); 
    dEnergies(ctr,2) = E_small(2)-E_small(3); 
    dEnergies(ctr,3) = E_small(3)-E_small(4); 
    dEnergies(ctr,4) = E_small(3)-E_small(5); 
    dEnergies(ctr,5) = E_small(2)-E_small(1); 
    plotQCL(Psi,E,Vx,x,0,1,0,num2str(b));
    axis([x(1)*1E10 x(end)*1E10 0 0.4]);
    Mov(ctr) = getframe;
    Energies(ctr,:) = E;
    ctr = ctr+1;

end

% %%
% indices = [7 6 5 4 3 2 1]
% % indices = [14 13 12 10 9 7 6] % if we simulate 3 or more modules. 
% 
% E_INJ_1= Energies(:,indices(1))*1E3; E_INJ_2= Energies(:,indices(2))*1E3;
% E_ULL= Energies(:,indices(3))*1E3; E_LLL_1= Energies(:,indices(4))*1E3;
% E_LLL_2= Energies(:,indices(5))*1E3; E_DEPOP_1 = Energies(:,indices(6))*1E3;
% E_DEPOP_2 = Energies(:,indices(7))*1E3;
% figure;
% % moduleBias = bias*54.8*1E3*1E-7;
% moduleBias = bias;
% hold on; 
% 
% plot(moduleBias,E_INJ_1-E_INJ_2,'-r','Linewidth',2.0);
% plot(moduleBias,E_INJ_1-E_ULL,'-g','Linewidth',2.0);
% plot(moduleBias,E_INJ_2-E_ULL,'-b','Linewidth',2.0);
% plot(moduleBias,E_ULL-E_LLL_1,'-y','Linewidth',2.0);
% plot(moduleBias,E_ULL-E_LLL_2,'-k','Linewidth',2.0);
% % ylim([0,6])
% 
% legend('\DeltaE_{INJ_1 INJ_2}','\DeltaE_{INJ_1 ULL }','\DeltaE_{INJ_2 ULL}','\DeltaE_{ULL LLL_1}', '\DeltaE_{ULL LLL_2}');
% 
% figure; hold on;
% plot(moduleBias,E_LLL_1-E_DEPOP_1,'-m','Linewidth',2.0);
% plot(moduleBias,E_LLL_1-E_DEPOP_2,'-k','Linewidth',2.0);
% 
% 
% legend('\DeltaE_{LLL_1 DEPOP_1}','\DeltaE_{LLL_1 DEPOP_2}');
% ylabel('\DeltaE [meV]'); xlabel('Bias kV/cm');
% 
% figure; hold on ; 
% plot(moduleBias,E_INJ_1, moduleBias,E_INJ_2,moduleBias,E_ULL); 
% 
% legend('E_{INJ1}','E_{INJ2}','E_{ULL}');
% ylabel('E [meV]'); xlabel('Bias kV/cm');
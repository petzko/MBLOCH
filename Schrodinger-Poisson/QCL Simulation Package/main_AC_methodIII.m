% This progrma will plot the conduction band edge potential and the solved
% moduli-squared wavefunctions for a quantum cascad lasers (QCLs). The QCL
% strucute is in a file in the folder "Input." Solved wavefunctions, energy
% levels, and other quantitites are saved in the folder "Output."
% This QCL simulation tool has been written by 
% Muhammad A. Talukder, Dept. of CSEE, UMBC
% email: anisuzzaman@umbc.edu
%--------------------------------------------------------------------------


close all;
bias = 0;
AC_energies_left = zeros(length(bias),2); 
AC_energies_right = zeros(length(bias),2); 
const = Constants();
nper = 2;
dEnergies = zeros(length(bias),4); 

ctr =1; 
for b = bias
        b
        efield.direction = '->';            % electric field direction
        efield.value = b;                  % electric field value
        indir = 'Input';
        infile = strcat(indir,'/','WELL.m');  % input file
        outdir = 'Output';
        nlevel = 2;                        % number of states to be solved
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

        const = Constants(); 
        %x in m and Vx in eV, bias in kV/cm 
        VB = Vx(1);
        
        if efield.direction == '->'
            Vtot = Vx-x*b*1E3*1E2;
        else
            Vtot = Vx+x*b*1E3*1E2;
        end
        
        % TB Part
        efield.value = b;                  % electric field value
        indir = 'Input';

        infile = strcat(indir,'/','WELLTB.m');  % input file
        outdir = 'Output';
        nlevel = 2                        % number of states to be solved
        T = 50;                            % Temperature
        nonpar_flag = 0;                    % non-parabolicity flag

        %--------------------------------------------------------------------------
        [ structure_TB,dipole_TB,FD_TB] = sbs(outdir,efield,infile,nlevel,T,nonpar_flag);
        TB_W = structure_TB(1,1)*1E-10; 
        %TB Basis confinement barrier width in Angstrom CB-TB barrier width difference; 
        dW = CB_W-TB_W; dL =(0*(nper-2)*Lp+dW); 
    
        
        ATB = load(strcat(outdir,'/','Psiout.txt'));
        E = load(strcat(outdir,'/','EnergyValues.txt'));
        V0 = load(strcat(outdir,'/','Vxmesh.txt'));
        dxTB = load(strcat(outdir,'/','GridSize.txt'));
        xTB = load(strcat(outdir,'/','x.txt'));
        nxTB = length(V0);
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
        
        
        numpts = round(Lp/dx); 
        for kk = nlevel+1 : 2*nlevel
            baselvl = PsiTB(:,kk-nlevel);
            shiftlvl = interp1(x+Lp,baselvl,x,'linear',0);
            PsiTB(:,kk) = shiftlvl; 
            ETB(kk) = ETB(kk-nlevel)+b*1E3*1E2*Lp;
        end
        
        VB = Vx(1); 
        Vtot = 0*Vx; 
        
        if efield.direction == '->'
            Vtot = Vx-x*b*1E3*1E2;
        else
            Vtot = Vx+x*b*1E3*1E2;
        end
        V0 = Vtot.*(x>(Lp+CB_W/2)).*(x<(2*Lp+CB_W/2));
        Vd = Vtot(1)*(x<=(Lp+CB_W/2))+ Vtot(1)*(x>=(2*Lp+CB_W/2));
       
        VxTB = x*b*1E3*1E2+V0+Vd;
        dV = Vtot-V0-Vd;

%         
% %         plot(x,dV);
%         
        [Hext,D,S] = calculateHamiltonian(+dV,x,PsiTB,ETB);
        Hnew = inv(D)*Hext; 
        [a,d] = eigs(Hnew,length(ETB),'sr'); 
        d = diag(d); 
%         [d,idx] = sort(d); 
%         a = a(:,idx);
        Psi_new = 0*PsiTB; 
        for j = 1:length(ETB)
            for i = 1:length(ETB)
                Psi_new(:,j) = Psi_new(:,j) + a(i,j)*PsiTB(:,i);
            end
            Psi_new(:,j) = Psi_new(:,j)/(sqrt(trapz(x,Psi_new(:,j).^2)));
        end
          
%       simply shift the wavefunctions from base to home with corresponding
%       energy and corresponding length! 
%       base = [2 6 7]; home_left = [1 3 4]; home_right = [5 8 9]
%       for kk = 1 : length(base)
%             baselvl = Psi_new(:,base(kk));
%             shiftlvl = interp1(x-Lp,baselvl,x,'linear',0);
%             Psi_new(:,home_left(kk)) = shiftlvl; 
%             d(home_left(kk)) = d(base(kk))-b*1E3*1E2*Lp;
%             
%             shiftlvl = interp1(x+Lp,baselvl,x,'linear',0);
%             Psi_new(:,home_right(kk)) = shiftlvl; 
%             d(home_right(kk)) = d(base(kk))+b*1E3*1E2*Lp;
%         end
%         
        plotQCL(Psi,E_ext,Vx,x,0,1,0,'ext');
        plotQCL(PsiTB,ETB,Vx,x,0,1,0,'TB');  
        plotQCL(Psi_new,d,Vx,x,0,1,0,'estimate');        
        
% 
end 



% This progrma will plot the conduction band edge potential and the solved
% moduli-squared wavefunctions for a quantum cascad lasers (QCLs). The QCL
% strucute is in a file in the folder "Input." Solved wavefunctions, energy
% levels, and other quantitites are saved in the folder "Output."
% This QCL simulation tool has been written by 
% Muhammad A. Talukder, Dept. of CSEE, UMBC
% email: anisuzzaman@umbc.edu
%--------------------------------------------------------------------------

clear;
close all;
bias = 11;
AC_energies = zeros(length(bias),2); 
const = Constants();
nper = 3;
dEnergies = zeros(length(bias),3); 
H_TB = zeros(7,7,length(bias)); 

ctr =1; 
for b = bias
        b
        efield.direction = '->';            % electric field direction
        efield.value = b;                  % electric field value
        indir = 'Input';
        infile = strcat(indir,'/','QCL183S.m');  % input file
        outdir = 'Output';
        nlevel = 9;                        % number of states to be solved
        T = 50;                            % Temperature
        nonpar_flag = 0;                    % non-parabolicity flag

        %--------------------------------------------------------------------------
        [structure_ext,dipole_ext,FD_ext ] =  sbs(outdir,efield,infile,nlevel,T,nonpar_flag);
        CB_W = structure_ext(1,1)*1E-10; %coupled basis confinement barrier width
        Lp = (sum(structure_ext(1,1:end-1))*1E-10)/nper; 
        
        
        A = load(strcat(outdir,'/','Psiout.txt'));
        E_ext = load(strcat(outdir,'/','EnergyValues.txt'));
        Egx = load(strcat(outdir,'/','Egxmesh.txt'));
        mx = load(strcat(outdir,'/','mxmesh.txt'));
        dx = load(strcat(outdir,'/','GridSize.txt'));
        Vx = load(strcat(outdir,'/','Vxmesh.txt'));
        x = load(strcat(outdir,'/','x.txt'));
        nx = length(Vx);

        Psi_ext = zeros(nx,nlevel);

        for kk = 1 : nlevel
            sindex = (kk-1)*nx + 1;
            findex = sindex+nx-1;
            Psi_ext(:,kk) = A(sindex : findex);
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
      
        
        % TB Part
        efield.value = b;                  % electric field value
        indir = 'Input';

        infile = strcat(indir,'/','QCL183STB.m');  % input file
        outdir = 'Output';
        nlevel = 5;                       % number of states to be solved
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
        Psi_old = zeros(nxTB,nlevel);
        Psi_periods = zeros(nx,nlevel,3);
        Psi_tot = zeros(nx,3*nlevel);
       
        E_periods = zeros(nlevel,3); 
        E_tot = zeros(3*nlevel,1); 
        
        
        % 1 <=> L 
        % 2 <=> 0 
        %   and 
        % 3 <=> R
        
        
        VpM = b*1E5*Lp; % voltage per module. 
        for kk = 1:nlevel
            sindex = (kk-1)*nxTB + 1;
            findex = sindex+nxTB-1;
            baselvl = ATB(sindex : findex);
            Psi_old(:,kk)=baselvl;
            baselvl = interp1(xTB+dL+Lp,baselvl,x,'linear',0);
            Psi_periods(:,kk,2)=baselvl;
            E_periods(kk,2) = E(kk)+b*1E3*1E2*(dL+Lp);
            
            %shift the remaining levles
            Psi_periods(:,kk,1) = interp1(x-Lp,Psi_periods(:,kk,2),x,'linear',0);
            E_periods(kk,1) = E_periods(kk,2)-b*1E3*1E2*Lp;
            
            Psi_periods(:,kk,3) = interp1(x+Lp,Psi_periods(:,kk,2),x,'linear',0);
            E_periods(kk,3) = E_periods(kk,2)+b*1E3*1E2*Lp;
            
            Psi_tot(:,kk) = Psi_periods(:,kk,1); E_tot(kk) = E_periods(kk,1); 
            Psi_tot(:,kk+nlevel) = Psi_periods(:,kk,2); E_tot(kk+nlevel) = E_periods(kk,2); 
            Psi_tot(:,kk+2*nlevel) = Psi_periods(:,kk,3); E_tot(kk+2*nlevel) = E_periods(kk,3); 
        end

        %x in m and Vx in eV, bias in kV/cm 
        VB = Vx(1);
        
        if efield.direction == '->'
            Vtot = Vx-x*b*1E3*1E2;
        else
            Vtot = Vx+x*b*1E3*1E2;
        end
        
        VL = Vtot.*(x>=(0)).*(x<(Lp+CB_W/2));
        V0 = Vtot.*(x>=(Lp+CB_W/2)).*(x<(2*Lp+CB_W/2));
        VR = Vtot.*(x>=(2*Lp+CB_W/2)).*(x<(3*Lp+CB_W));
        Vd_L = Vtot(1)*(x>=Lp+CB_W/2);
        Vd_0 = Vtot(1)*(x<(Lp+CB_W/2))+ Vtot(1)*(x>=(2*Lp+CB_W/2));
        Vd_R = Vtot(1)*(x<(2*Lp+CB_W/2));
       
        dV(:,1) = VR+V0-Vd_L;
        dV(:,2) = VL+VR-Vd_0;
        dV(:,3) = VL+V0-Vd_R;
        
        [Hext,D,S] = calculateHamiltonianNEW(Psi_periods,E_periods,dV,x);
        VTB = Vtot(1)*(x<(Lp+CB_W))+ Vtot(1)*(x>=(2*Lp))+x*1E5*bias;
        
        indices01 = [7 6 5 4 3 2 1];
        indices02 = [12 11 10 9 8 7 6];
        
        for i = 1:length(indices01)
            idx_i = indices01(i);
            for j = 1:length(indices01)
                idx_j = indices01(j);
                if( i==j)
                    H_TB(i,j,ctr) = Hext(idx_i,idx_j);
                else
                    H_TB(i,j,ctr) = S(idx_i,idx_j); 
                end
             end
        end
        
        
        INJ1 = 7 ; INJ2 = 6; ULL = 5; LLL1 = 4; LLL2 =3;
        
        dEnergies(ctr,1) = Hext(INJ1,INJ1)-Hext(ULL,ULL); 
        dEnergies(ctr,2) = Hext(INJ2,INJ2)-Hext(ULL,ULL); 
        dEnergies(ctr,3) = Hext(INJ1,INJ1)-Hext(INJ2,INJ2); 
        dEnergies(ctr,4) = Hext(ULL,ULL)-Hext(LLL1,LLL1); 
        dEnergies(ctr,5) = Hext(ULL,ULL)-Hext(LLL2,LLL2); 

        
              
       
        AC_energies(ctr,1) = S(INJ1,ULL); 
%         AC_energies(ctr,2) = S(ULL,INJ1); 
        AC_energies(ctr,3) = S(INJ2,ULL); 
%         AC_energies(ctr,4) = S(ULL,INJ2); 
        
        
        
%        Hnew = inv(D)*Hext; 
        [a,d] = eigs(Hext,D,length(Hext),'sr'); 
        d = diag(d); 
        
            
        Psi_new = 0*Psi_tot; 
        for j = 1:3*nlevel
            for i = 1:nlevel
                Psi_new(:,j) = Psi_new(:,j) + a(i,j)*Psi_periods(:,i,1) +a(i+nlevel,j)*Psi_periods(:,i,2) +a(i+2*nlevel,j)*Psi_periods(:,i,3) ;
            end
            Psi_new(:,j) = Psi_new(:,j)/(sqrt(trapz(Psi_new(:,j).^2)));
        end
        
 
          ctr = ctr + 1; 
        plotQCL(Psi_tot,E_tot,[Vx,VTB],x,0,1,0,'TB');  
%         plotQCL(Psi_ext,E_ext,Vx,x,0,1,0,'Coupled');  
%         plotQCL(Psi_new,d,Vx,x,0,1,0,'TB->Coupled');  
        display(['HTB @ bias: ' num2str(b)]); 
        display(num2str(reshape(H_TB',1,49)));

end 

val_field = b*1e3*1e2;
V_TB = val_field*x + Vd_0;


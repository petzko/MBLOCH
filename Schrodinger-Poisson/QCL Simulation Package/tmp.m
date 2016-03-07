   % TB Part

        efield.value = b;                  % electric field value
        efield.direction = '->';            % electric field direction
        indir = 'Input';

        infile = strcat(indir,'/','QCL183STB.m');  % input file
        outdir = 'Output';
        nlevel = 5;                        % number of states to be solved
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

        PsiTB = zeros(nx,3*nlevel);
        ETB = zeros(3*nlevel,1); 
        %length of one period. 
        for kk = 1:nlevel
            sindex = (kk-1)*nx + 1;
            findex = sindex+nx-1;
            PsiTB(:,kk) = A(sindex : findex);
            ETB(kk) = E(kk);
        end
        % plotQCL(PsiTB,ETB,Vx,x,0,1,0,'TB'); 

        % right WFS
        Lp = 548E-10; VpM = b*1E5*Lp; % voltage per module. 

        numpts = round(Lp/dx); 
        for kk = nlevel+1 : 2*nlevel
            baselvl = PsiTB(:,kk-nlevel);
            shiftlvl = interp1(x+Lp,baselvl,x,'linear',0);
            PsiTB(:,kk) = shiftlvl; %shiftlvl;
            ETB(kk) = E(kk-nlevel)+b*1E3*1E2*Lp;
        end

        %left WFS

        for kk = 2*nlevel+1 : 3*nlevel
            baselvl = PsiTB(:,kk-2*nlevel);
            shiftlvl = interp1(x-Lp,baselvl,x,'linear',0);
            PsiTB(:,kk) = shiftlvl;
            ETB(kk) = E(kk-2*nlevel)-b*1E3*1E2*Lp;
        end

        plotQCL(PsiTB(:,1:nlevel),ETB(1:nlevel),Vx,x,0,1,0,'mid');
%         plotQCL(PsiTB(:,nlevel+1:2*nlevel),ETB(nlevel+1:2*nlevel),Vx,x,0,1,0,'right');
%         plotQCL(PsiTB(:,2*nlevel+1:3*nlevel),ETB(2*nlevel+1:3*nlevel),Vx,x,0,1,0,'left');
%         plotQCL(PsiTB,ETB,Vx,x,0,1,0,'total');
        hold on; 
%         plot(x*1E10,VB-V0);
        dV = VB-V0; 

        INJ2 = 1; ULL = 5; 

        AC_right = zeros(nlevel); 
        AC_left = zeros(nlevel); 
        for i = 1:nlevel
            for j =  1:nlevel
                Psi_i = PsiTB(:,i); 
                Psi_j = PsiTB(:,j+nlevel);   
                %do not multiply by dx since the Psi's are normalized such that
                % trapz(Psi_i^2) = 1
                AC_right(i,j) = trapz(conj(Psi_i).*dV.*conj(Psi_j));
            end
        end
        dE1 = AC_right(ULL,INJ2);dE2 = AC_left(INJ2,ULL); 
        for i = 1:nlevel
            for j =  1:nlevel
                Psi_i = PsiTB(:,i); 
                Psi_j = PsiTB(:,j+2*nlevel);        
                AC_left(i,j) = trapz(conj(Psi_i).*dV.*conj(Psi_j));
            end
        end
        
        AC_energies(ctr) = abs(dE1);
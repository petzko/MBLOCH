function [x,nx,Egx,Vx,mx,dopx,rindexx,epsr,hbarwlo,Vxflat] = QCLMesh(outdir,structure,dx,efield,T,purge)

% This function creates the mesh for bandgap, potential,
% electron effective mass, doping density, and effective index.
%
% USAGE:
%
% [x,nx,Egx,Vx,mx,dopx,rindexx] = QCLMesh(structure,dx,efield)
%
% INPUT:
%
% structure - the heterostructure superlattice
%             1st column of 'structure' is the thicknesses of the layers
%             2nd column of 'structure' is the material
% dx - grid spacing
% efield - externally applied electric field
%   efield.value - value of the externally applied electric field
%   efield.direction - direction of the efield
% T - temperature
%
% OUTPUT:
%
% x - vector specifying mesh coordinates
% nx - size of the potential mesh
% Egx - bandgap potential mesh ( Eg(x) )
% Vx - potential mesh ( V(x) )
% mx - mass mesh ( m(x) )
% dopx - doping density mesh
% rindexx - refractive index mesh
%
% Author - Muhammad Anisuzzaman Talukder (anisuzzaman@umbc.edu)
%--------------------------------------------------------------

const = Constants();

% Following section of the code selects bandgap,
% conduction band offset, effective mass, and
% doping density according to the material system.
% Then meshs for each of the above quantities are
% generated for the quantum structure.
% 'structure' input variable has four columns:
%   Column 1: thickness of the layers
%   Column 2: material system
%       '0' is for GaAs
%       '1' is for GaAlAs
%       '2' is for InGaAs
%       '3' is for InAlAs
%   Column 3: y - composition parameter (amount of Al/Ga in InAs or GaAs)
%   Column 4: doping density

structure(1,:) = structure(1,:) * 1e-10;
nlayers = length(structure(1,:));
nx = 0;
sindex = 1;

%%%% params PeTz

cob=0.135;
EgGaAs=1.519-0.5405e-3*T^2/(T+204); EgInAs=0.417-0.276e-3*T^2/(T+94);
EgAlAs=3.099-0.885e-3*T^2/(T+530); DsoGaAs=0.341;     
DsoInAs=0.39; DsoAlAs=0.28;
FGaAs=-1.94; FInAs=-2.9;
FAlAs=-0.48; EpGaAs=28.8;
EpInAs=21.5; EpAlAs=21.1;
%%%% end params PeTz

for ii = 1 : nlayers
    
    nxlayer(ii) = round(structure(1,ii)/dx);
    
    material = structure(2,ii);
    y = structure(3,ii);
    
    switch material
 %%%%%%%%%%%%%%%%%% begin PetZ %%%%%%%%
% %   case 0
% %             Eg = 1.519 - 5.41e-4*T*T/(T + 204); % Ref. 1
% %             V = 0;
% %             m = 0.063 * const.me; % Ref. 1
% %             rindex = 3.3; % Ref. 1 
% %             epsr.high = 10.89; % Ref. 1
% %             epsr.static = 12.9; % Ref. 1
% %             hbarwlo = 0.036; % Ref. 1
% %             
% %         case 1
% %             Eg = 1.519 + 1.155*y + 0.37*y*y - 5.41e-4*T*T/(T + 204); % Ref. 1
% %             V = 0.62 * (1.155*y + 0.37*y*y);
% %             m = (0.063 + 0.083*y) * const.me; % Ref. 1
% %             rindex = 3.3 - 0.53*y + 0.09*y*y; % Ref. 1
            
        case 0
            
           V=0;
           m=1/(1+2*FGaAs+EpGaAs*(EgGaAs+2*DsoGaAs/3)/EgGaAs/(EgGaAs+DsoGaAs))*const.me; 
           Eg=EgGaAs;
           rindex = 3.3; % Ref. 1 
           epsr.high = 10.89; % Ref. 1
           epsr.static = 12.9; % Ref. 1
           hbarwlo = 0.036; % Ref. 1

        case 1
            V=cob;
            m=y/(1+2*FAlAs+EpAlAs*(EgAlAs+2*DsoAlAs/3)/EgAlAs/(EgAlAs+DsoAlAs));
            m=m+(1-y)/(1+2*FGaAs+EpGaAs*(EgGaAs+2*DsoGaAs/3)/EgGaAs/(EgGaAs+DsoGaAs));
            m = m*const.me;
            Eg=EgGaAs*(1-y)+EgAlAs*y+y*(1-y)*0.127 ;
            rindex = 3.3 - 0.53*y + 0.09*y*y; % Ref. 1
%%%%%%%%%%%%%%%%%% end PetZ %%%%%%%%


        case 2
            Eg = 0.42 + 0.625*y - (5.8/(T+300)-4.19/(T+271)) * 1e-4*T^2*y -...
                 4.19*1e-4*(T^2/(T+271)) + 0.475*y^2; % Ref. 1
            %Eg = 0.816 - 4e-4*T*T/(T + 145);
            V = 0;
            %----------------------
            if y == 0.1 && structure(2,ii-1) == 2
                V = 0.52-0.7526;
            end
            %----------------------
            m = (0.023 + 0.037*y + 0.003*y^2) * const.me; % Ref. 1
            rindex = 3.51 - 0.16*y; % Ref. 1
            epsr.high = 12.3 - 1.4*y; % Ref. 1
            epsr.static = 15.1 - 2.87*y + 0.67*y^2; % Ref. 1
            hbarwlo = 0.034; % this value is for Ga(0.47)Al(0.53)As          
            
        case 3
            %Eg = 0.816 + 0.65 - 4e-4*T*T/(T + 145);
            Eg = 0.36 + 2.012*y + 0.698*y^2; % Ref. 2
            
            if ii == 1
                 if structure(2,ii+1) == 2
                    yi = structure(3,ii+1);
                else
                    if structure(2,ii+2) == 2
                        yi = structure(3,ii+2);
                    else
                        yi = structure(3,ii+3);
                    end
                end
                
            else
                if structure(2,ii-1) == 2
                    yi = structure(3,ii-1);
                else
                    if structure(2,ii+1) == 2
                        yi = structure(3,ii+1);
                    else
                        yi = structure(3,ii+2);
                    end
                end
            end

            Egi = 0.42 + 0.625*yi - (5.8/(T+300)-4.19/(T+271)) * 1e-4*T^2*yi -...
                 4.19*1e-4*(T^2/(T+271)) + 0.475*yi^2; % Eg of InGaAs
            V = (0.653 + 0.1*(1-y))*(Eg-Egi); % Ref. 1
            m = (0.0427 + 0.0333) * const.me;            
            rindex = 3.51 - 0.16*y - 0.287;
    end % end switch material
        
    findex = sindex + nxlayer(ii) - 1;
    
    Egx(sindex : findex) = Eg;
    Vx(sindex : findex) = V;
    mx(sindex : findex) = m;
    dopx(sindex : findex) = structure(4,ii)*1e23;
    rindexx(sindex : findex) = rindex;
    
    sindex = sindex + nxlayer(ii);
    
    nx = nx + nxlayer(ii);
end

if structure(2,1) == 2 || structure(2,1) == 3
    Vx = CalcBandOffset(structure,dx);
end

if purge.flag == 1
    A1 = max(Vx);
    A2 = min(Vx);
    B1 = max(Egx);
    B2 = min(Egx);
    C1 = max(mx);
    C2 = min(mx);
    D1 = min(rindexx);
    D2 = max(rindexx);
    
    Vx = [];
    Vx = zeros(1,nx);
    Egx = [];
    Egx = zeros(1,nx);
    mx = [];
    mx = zeros(1,nx);
    rindexx = [];
    rindexx = zeros(1,nx);
    
    sindex = 1;
   for ii = 1 : nlayers
       layerVx = [];
       layerEgx = [];
       layermx = [];
       layerrindexx = [];
       material = structure(2,ii);

       switch material
           case 0
               V1 = A1;
               V2 = A2;
               Eg1 = B1;
               Eg2 = B2;
               M1 = C1;
               M2 = C2;
               R1 = D1;
               R2 = D2;
               if purge.time(1) >= purge.maxtime
                   purge1 = 1e-12;
               else
                   purge1 = purge.npdepth*(1 - purge.time(1)/purge.maxtime);
               end
               if purge.time(2) >= purge.maxtime
                   purge2 = 1e-12;
               else
                   purge2 = purge.npdepth*(1 - purge.time(2)/purge.maxtime);
               end
           case 1
               V1 = A2;
               V2 = A1;
               Eg1 = B2;
               Eg2 = B1;
               M1 = C2;
               M2 = C1;
               R1 = D2;
               R2 = D1;
               if purge.time(1) >= purge.maxtime
                   purge1 = 1e-12;
               else
                   purge1 = purge.npdepth*(1 - purge.time(2)/purge.maxtime);
               end
               if purge.time(2) >= purge.maxtime
                   purge2 = 1e-12;
               else
                   purge2 = purge.npdepth*(1 - purge.time(1)/purge.maxtime);
               end
           case 2
               V1 = A1;
               V2 = A2;
               Eg1 = B1;
               Eg2 = B2;
               M1 = C1;
               M2 = C2;
               R1 = D1;
               R2 = D2;
               if purge.time(1) >= purge.maxtime
                   purge1 = 1e-12;
               else
                   purge1 = purge.npdepth*(1 - purge.time(1)/purge.maxtime);
               end
               if purge.time(2) >= purge.maxtime
                   purge2 = 1e-12;
               else
                   purge2 = purge.npdepth*(1 - purge.time(2)/purge.maxtime);
               end
           case 3
               V1 = A2;
               V2 = A1;
               Eg1 = B2;
               Eg2 = B1;
               M1 = C2;
               M2 = C1;
               R1 = D2;
               R2 = D1;
               if purge.time(1) >= purge.maxtime
                   purge1 = 1e-12;
               else
                   purge1 = purge.npdepth*(1 - purge.time(2)/purge.maxtime);
               end
               if purge.time(2) >= purge.maxtime
                   purge2 = 1e-12;
               else
                   purge2 = purge.npdepth*(1 - purge.time(1)/purge.maxtime);
               end
       end % end switch material

       for kk = 1 : nxlayer(ii)
         %  if kk <= nxlayer(ii)/2
         %      if ii == 1
         %          layerVx(kk) = V2;
         %          layerEgx(kk) = Eg2;
         %          layermx(kk) = M2;
         %          layerrindexx(kk) = R2;
                   
         %      else
         %          layerVx(kk) = (V1+V2)/2 + (V2-V1)/2*tanh(((kk-1)*dx-0)/purge1);
         %          layerEgx(kk) = (Eg1+Eg2)/2 + (Eg2-Eg1)/2*tanh(((kk-1)*dx-0)/purge1);
         %          layermx(kk) = (M1+M2)/2 + (M2-M1)/2*tanh(((kk-1)*dx-0)/purge1);
         %          layerrindexx(kk) = (R1+R2)/2 + (R2-R1)/2*tanh(((kk-1)*dx-0)/purge1);
         %      end
         %  else
               if ii == nlayers
                   layerVx(kk) = V2;
                   layerEgx(kk) = Eg2;
                   layermx(kk) = M2;
                   layerrindexx(kk) = R2;
               else
                   layerVx(kk) = (V1+V2)/2 + (V1-V2)/2*tanh(((kk-1)*dx-structure(1,ii))/purge2);
                   layerEgx(kk) = (Eg1+Eg2)/2 + (Eg1-Eg2)/2*tanh(((kk-1)*dx-structure(1,ii))/purge2);
                   layermx(kk) = (M1+M2)/2 + (M1-M2)/2*tanh(((kk-1)*dx-structure(1,ii))/purge2);
                   layerrindexx(kk) = (R1+R2)/2 + (R1-R2)/2*tanh(((kk-1)*dx-structure(1,ii))/purge2);
               end
          % end
       end
       findex = sindex + nxlayer(ii) - 1;
       Vx(sindex:findex) = layerVx;    
       Egx(sindex:findex) = layerEgx;
       mx(sindex:findex) = layermx;
       rindexx(sindex:findex) = layerrindexx;
       sindex = sindex + nxlayer(ii);
   end % end for ii = 1 : nlayers
end % end if purge_flag == 1

VCB = max(Vx)-min(Vx);
x = (0 : (nx-1))'*dx;

efield.value = efield.value * 1e3 * 1e2;

switch efield.direction
    case '->'
        vx = (1:1:nx)*dx*efield.value;
    case '<-'
        vx = (nx:-1:1)*dx*efield.value;
    otherwise
        error('Direction of the applied electric field is not set properly'); 
end

Vxflat = Vx;    % potential mesh without field
Vx = Vx + vx;   % potential mesh with field

% Writing the parameter values to the output file 'Output.txt'
fid = fopen(strcat(outdir,'/','SolveBandStructure.txt'),'a');
fprintf(fid,'\n%s\n','Material parameters');
fprintf(fid,'%s\n','-------------------');
if nlayers > 1
    noiter = 2;
else
    noiter = 1;
end

for jj = 1 : noiter
    
    material = structure(2,jj);
    y = structure(3,jj);
    
    switch material
 
        case 0
            Eg = 1.519 - 5.41e-4*T*T/(T + 204);
            V = 0;
            m = 0.067 * const.me;
            rindex = 3.3;
            if jj == 1
                matindex = 1;
            else
                matindex = 2;
            end
            fprintf(fid,'%s %1.0f %s\n','Material',matindex,': GaAs');
            fprintf(fid,'    %s %1.2f %s %1.2f\n','Ga:',1-y,'Al',y);
            fprintf(fid,'    %s %2.2e\n','m*(GaAs):',m);
            fprintf(fid,'    %s %2.4f %s\n','Eg(GaAs):',Eg,'eV'); 
            fprintf(fid,'    %s %2.4f\n','Refractive index (GaAs):',rindex);
            fprintf(fid,'    %s\n','Dielectric Constant');
            fprintf(fid,'      %s %2.2f\n','static:',epsr.static);
            fprintf(fid,'      %s %2.2f\n','high frequency:',epsr.high);
            fprintf(fid,'    %s %2.2f\n\n','Optical phonon energy (meV):',1e3*hbarwlo);
        case 1
            Eg = 1.519 + 1.155*y + 0.37*y*y - 5.41e-4*T*T/(T + 204);
            V = 0.62 * (1.155*y + 0.37*y*y);
            m = (0.067 + 0.083*y) * const.me;
            rindex = 3.3 - 0.53*y + 0.09*y*y;
            if jj == 1
                matindex = 1;
            else
                matindex = 2;
            end
            fprintf(fid,'%s %1.0f %s\n','Material',matindex,': AlGaAs');
            fprintf(fid,'    %s %1.2f %s %.1.2f\n','Ga:',1-y,'Al:',y);
            fprintf(fid,'    %s %2.2e\n','m*(AlGaAs):',m);
            fprintf(fid,'    %s %2.4f %s\n','Eg(AlGaAs):',Eg,'eV'); 
            fprintf(fid,'    %s %2.4f\n\n','Refractive index (AlGaAs):',rindex);
                        
        case 2
            Eg = 0.42 + 0.625*y - (5.8/(T+300)-4.19/(T+271)) * 1e-4*T^2*y -...
                 4.19*1e-4*(T^2/(T+271)) + 0.475*y^2; % Ref. 1
            %Eg = 0.816 - 4e-4*T*T/(T + 145);
            V = 0;
            m = (0.023 + 0.037*y + 0.003*y^2) * const.me; % Ref. 1
            rindex = 3.51 - 0.16*y; % Ref. 1
            if jj == 1
                matindex = 1;
            else
                matindex = 2;
            end
            fprintf(fid,'%s %1.0f %s\n','Material',matindex,': InGaAs');
            fprintf(fid,'    %s %1.2f %s %1.2f\n','In:',1-y,'Ga:',y);
            fprintf(fid,'    %s %2.2e\n','m*(InGaAs):',m);
            fprintf(fid,'    %s %2.4f %s\n','Eg(InGaAs):',Eg,'eV'); 
            fprintf(fid,'    %s %2.4f\n','Refractive index (InGaAs):',rindex);
            fprintf(fid,'    %s\n','Dielectric Constant');
            fprintf(fid,'      %s %2.4f\n','static:',epsr.static);
            fprintf(fid,'      %s %2.4f\n','high frequency:',epsr.high);
            fprintf(fid,'    %s %2.2f\n\n','Optical phonon energy (meV):',1e3*hbarwlo);
                        
        case 3
            %Eg = 0.816 + 0.65 - 4e-4*T*T/(T + 145);
            Eg = 0.36 + 2.012*y + 0.698*y^2; % Ref. 2
            
            if ii == 1
                yi = structure(3,ii+1);
            else
                yi = structure(3,ii-1);
            end

            Egi = 0.42 + 0.625*yi - (5.8/(T+300)-4.19/(T+271)) * 1e-4*T^2*yi -...
                 4.19*1e-4*(T^2/(T+271)) + 0.475*yi^2; % Eg of InGaAs
            V = (0.653 + 0.1*(1-y))*(Eg-Egi); % Ref. 1
            m = (0.0427 + 0.0333) * const.me;            
            rindex = 3.51 - 0.16*y - 0.287;
            if jj == 1
                matindex = 1;
            else
                matindex = 2;
            end
            fprintf(fid,'%s %1.0f %s\n','Material',matindex,': InAlAs');
            fprintf(fid,'    %s %1.2f %s %1.2f\n','In:',1-y,'Al:',y);
            fprintf(fid,'    %s %2.2e\n','m*(InAlAs):',m);
            fprintf(fid,'    %s %2.4f %s\n','Eg(InAlAs):',Eg,'eV'); 
            fprintf(fid,'    %s %2.4f\n\n','Refractive index (InAlAs):',rindex);
            
                        
    end % end switch material
        
    
end

fprintf(fid,'%s %-2.4f %s\n','Conduction band offset',VCB,'eV');
fprintf(fid,'\n%s\n','Simulation parameters');
fprintf(fid,'%s\n','---------------------');
fprintf(fid,'%s %-6.0f\n','Number of points:',nx);
fprintf(fid,'%s %-2.2f %s\n','Grid size:',1e10*dx,'Angstrom');

fclose(fid);


% Reference
% [1] M. Levinshtein, S. Rumyantsev, and M. Shur, "Handbood Series on 
%     Semiconductor Parameters," World Scientific Publishing Co.,
%     ch. 3, 1999.
% [2] V. Mitin, V. Kochelap, and M. Stroscio, " Quantum Heterostructures,"
%     Cambridge Univ. Press, p-133, 1999.
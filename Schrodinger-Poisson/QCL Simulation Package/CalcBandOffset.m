function V = CalcBandOffset(structure,dx)
%
% This function calculates the band offset at InGaAs/InAlAs material
% interface
%
% USAGE:
%
% [V] = CalcStrain(infile,dx)
%
% INPUT:
%
% infile - the name of the input file that contains the QCL structure
% dx - grid size
%
% OUTPUT:
%
% V - the band offsets
%
% This program has been written on: July 11, 2009
%
% Author - Muhammad Anisuzzaman Talukder (anisuzzaman@umbc.edu)
%--------------------------------------------------------------------------
%fid = fopen(infile, 'r');
%structure = fscanf(fid, '%g ', [4 inf]);    
%fclose(fid);
%structure(1,:) = structure(1,:) * 1e-10;
nlayers = length(structure(1,:));
sindex = 1;

a_InAs = 6.0583;
a_GaAs = 5.6533;
a_AlAs = 5.6611;

G_InAs = 1.587;
G_GaAs = 2.522;
G_AlAs = 2.656;

D_InAs = 1.088;
D_GaAs = 0.934;
D_AlAs = 0.854;

Ec_InAs = -6.13;
Ec_GaAs = -5.29;
Ec_AlAs = -4.27;

ac_InAs = -5.08;
ac_GaAs = -7.17;
ac_AlAs = -5.64;

a_inplane = zeros(1,nlayers-1);
material = zeros(1,2);
y = zeros(1,2);
a = zeros(1,2);
G = zeros(1,2);
D = zeros(1,2);
eps_x = zeros(1,2);
eps_y = zeros(1,2);
eps_z = zeros(1,2);
h = zeros(1,2);
a_xplane = zeros(1,2);
Ec = zeros(1,2);
ac = zeros(1,2);
CF = zeros(1,2);

for ii = 1 : nlayers-1
    nxlayer = round(structure(1,ii)/dx);
    for jj = 1 : 2
        material(jj) = structure(2,ii+jj-1);
        y(jj) = structure(3,ii+jj-1);
        
        switch material(jj)
            case 2
                h(jj) = 6;
                a(jj) = (1-y(jj))*a_InAs + y(jj)*a_GaAs;
                G(jj) = (1-y(jj))*G_InAs + y(jj)*G_GaAs;
                D(jj) = (1-y(jj))*D_InAs + y(jj)*D_GaAs;
                Ec(jj) = (1-y(jj))*Ec_InAs + y(jj)*Ec_GaAs + ...
                         3*(1-y(jj))*y(jj)*(-ac_InAs+ac_GaAs)* ...
                         (a_InAs-a_GaAs)/a(jj);
                ac(jj) = (1-y(jj))*ac_InAs + y(jj)*ac_GaAs;
                CF(jj) = (0.47-y(jj)) * 0.529 * 0.1275;
            
            case 3
                h(jj) = 4;
                a(jj) = (1-y(jj))*a_InAs + y(jj)*a_AlAs;
                G(jj) = (1-y(jj))*G_InAs + y(jj)*G_AlAs;
                D(jj) = (1-y(jj))*D_InAs + y(jj)*D_AlAs;
                Ec(jj) = (1-y(jj))*Ec_InAs + y(jj)*Ec_AlAs + ...
                         3*(1-y(jj))*y(jj)*(-ac_InAs+ac_AlAs)* ...
                         (a_InAs-a_AlAs)/a(jj);
                ac(jj) = (1-y(jj))*ac_InAs + y(jj)*ac_AlAs;
                CF(jj) = (y(jj)-0.48) * 0.427 * 0.1275;
        end
         
    end
    
    a_inplane(ii) = (a(1)*G(1)*h(1)+a(2)*G(2)*h(2))/(G(1)*h(1)+G(2)*h(2));
    
    for jj = 1 : 2
        eps_x(jj) = a_inplane(ii)/a(jj) - 1;
        eps_y(jj) = eps_x(jj);
        a_xplane(jj) = a(jj)*(1-D(jj)*(a_inplane(ii)/a(jj)-1));
        eps_z(jj) = a_xplane(jj)/a(jj) - 1;
        Tr_eps = eps_x(jj)+eps_y(jj)+eps_z(jj);
        Ec(jj) = Ec(jj) + ac(jj)*Tr_eps;
    end
    
    dEc=Ec(2)-Ec(1);
    CF = 0.057 + CF(1) + CF(2);
    
    findex = sindex + nxlayer - 1;
    if structure(2,ii) == 2
        V(sindex:findex) = 0;
    end
    if structure(2,ii) == 3
        V(sindex:findex) = abs(dEc) - CF;
    end
    
    if ii == nlayers-1
        sindex = sindex + nxlayer;
        nxlayer = round(structure(1,ii+1)/dx);
        findex = sindex + nxlayer - 1;
        if structure(2,ii+1) == 2
            V(sindex:findex) = 0;
        end
        if structure(2,ii+1) == 3
            V(sindex:findex) = abs(dEc) - CF;
        end
    end
    
    sindex = sindex + nxlayer;
    
end
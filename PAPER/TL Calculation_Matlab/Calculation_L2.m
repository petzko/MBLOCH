clc; close all; clear;

%***********Assuption: infinite current sheet (2D) ******************
width = 60e-6;     % unit: m in x direction
I = 1;             % unit: A
J = I/width;       % unit: A/m
Jr = J*[0 1 0];    % unit: A/m in y direction
GridsX = 30;       % number of grids in x direction
GridsZ = 30;

% Parameters
mu0 = 1.2566370614e-6; %H/m
mu_r = 1;
% Central position of device [0,0,0]
Xmin=-width/2;
Xmax=width/2;
Zmin=Xmin;
Zmax=Xmax;
PlotXmin=Xmin*4;
PlotXmax=Xmax*4;
PlotZmin=PlotXmin;
PlotZmax=PlotXmax;
dx=(Xmax-Xmin)/GridsX; %step in the x direction
dz=(Zmax-Zmin)/GridsZ;

Bx=zeros(GridsX,GridsZ); % x component of field 
Bz=zeros(GridsX, GridsZ);% z component of field
[XData,ZData]=meshgrid(linspace(PlotXmin,PlotXmax,GridsX), linspace(PlotZmin,PlotZmax,GridsZ)); %build arrays of plot plane 


for m=1:GridsX     % in the x direction
    for n=1:GridsZ % in the z direction 
        PlotX=XData(m,n); 
        PlotZ=ZData(m,n);
        % on the current sheet
        if ((PlotZ==0)&&(PlotX>=Xmin)&&(PlotX<=Xmax))
            Bx(m,n)=0.5*mu0*Jr(2)*dx;
            Bz(m,n)=0;
            continue;
        end
        Rp=[PlotX 0 PlotZ]; % observed point
        for i=1:GridsX      % repeat in the x direction
        XCellCenter=Xmin+(i-1)*dx+0.5*dx; % center of current subsection 
        Rc=[XCellCenter 0 0]; % position vector of center of current subsection 
        R=Rp-Rc; 
        norm_R=norm(R);
        R_Hat=R/norm_R;
        dB=mu0*cross(Jr*dx,R_Hat)/(2*pi*norm_R); % Biot-Savart law in 2D
        Bx(m,n)=Bx(m,n)+dB(1,1); % add the x component
        Bz(m,n)=Bz(m,n)+dB(1,3); % add the z component
        end
    end
end
% Calculation of distributed inductance L' using Psi=L'*I=trapz(B,l)
IntCurve_l = round(GridsX/6);   %left point (bottom point)
IntCurve_r = GridsX-IntCurve_l; %right point (upper point)
Psi=(trapz(abs(Bz(IntCurve_l:IntCurve_r,IntCurve_l)))+trapz(abs(Bz(IntCurve_l:IntCurve_r,IntCurve_r)))...
    +trapz(abs(Bx(IntCurve_l,IntCurve_l:IntCurve_r)))+trapz(abs(Bx(IntCurve_r,IntCurve_l:IntCurve_r))))*dx;

L = Psi/I;
display(['L= ' num2str(L) ' H/m']);

% Plot magnetic field of cross section
quiver(XData*1e6, ZData*1e6, Bx, Bz);
title('Magnetic field in cross section');
xlabel('y(um)');
ylabel('z(um)');


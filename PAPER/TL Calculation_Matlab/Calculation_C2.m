clc; close all; clear;


%---------------- Enter the dimensions and requied iterations ------------%
width = 60;   % cavity width in um
height = 10;  % cavity thickness in um
Ni = 300;     % Number of iterations for the Poisson solver
U = 1;        % unit: V


% Parameters
Eps0 = 8.85e-12;     % unit F/m
Eps_r = 12.96;       % for active region, for air is 1

Nx = width*3;        % Length of plate in terms of number of X-grids
Ny = height*15;      % Length of plate in terms of number of Y-grids
mpx = ceil(Nx/2);    % Mid-point of x
lp = floor(width/2); % Mid point of plate
V = zeros(Nx,Ny);    % Potential (Voltage) matrix

% set boundary condition
V(1,:) = 0;
V(Nx,:) = 0;
V(:,1) = 0;
V(:,Ny) = 0;


for z = 1:Ni    % Number of iterations
        
        for i=2:Nx-1
        for j=2:Ny-1
            
            % The next two lines are meant to force the matrix to hold the 
            % potential values for all iterations
              V(mpx-lp:mpx+lp,height) = 1;
              V(mpx-lp:mpx+lp,1) = 0;
            % Helmholtz equation  
              V(i,j)=0.25*(V(i+1,j)+V(i-1,j)+V(i,j+1)+V(i,j-1));
        end
        end
        
end

% Take transpose for proper x-y orientation
V = V';

[Ex,Ey]=gradient(V);
Ex = -Ex;
Ey = -Ey;

% Electric field Magnitude
E = sqrt(Ex.^2+Ey.^2);  

x = (1:Nx)-mpx;
y = 1:Ny;

% Calculation of distributed capacitance C' using C=Qs/U
EEx=abs(Ex');
EEy=abs(Ey');
EE=E';
Qs=(trapz(EEy((mpx-2*lp):(mpx+2*lp),round(height/2)))*12.96+trapz(EEy((mpx-2*lp):(mpx+2*lp),height*2))+...
    trapz(EEx(mpx-2*lp,(round(height/2):height*2)))+trapz(EEx(mpx+2*lp,(round(height/2):height*2))))*Eps0;
C = Qs/U;
display(['C= ' num2str(C) ' F/m']);

% Contour Display for electric potential
figure(1)
contour_range_V = 0:0.01:1;
contour(x,y,V,contour_range_V,'linewidth',0.5);
axis([min(x) max(x) min(y) max(y)]);
colorbar('location','eastoutside','fontsize',14);
xlabel('y-axis in um','fontsize',14);
ylabel('z-axis in um','fontsize',14);
title('Electric Potential distribution, V(y,z) in volts','fontsize',14);
h1=gca;
set(h1,'fontsize',14);
fh1 = figure(1); 
set(fh1, 'color', 'white')

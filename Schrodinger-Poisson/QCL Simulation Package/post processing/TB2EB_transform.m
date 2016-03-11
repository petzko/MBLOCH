
%% define the QCL structure
hbar= 1.054571726E-34; %J per sec
e0 = 1.60217657E-19; % electron charge..
%barrier height [eV]
hb = 0.135; 
%bias in kV/cm converted to V/nm where the positive sign means growing in positive x direction
% Biasbase = +10.2;
Bias = 1E3*Biasbase/1E7; 

% extended hamiltonian CB profile:
Bar = [2.4 3.3 4.0 2.3 2.4];
Well = [9.2 16.2 7.4 7.6];

% slab profile starting barrier well barrier well
slab = [Bar(1) Well(1) Bar(2) Well(2) Bar(3) Well(3) Bar(4) Well(4) Bar(5)]; 

Lp = 0;
for i=1:length(slab)
    Lp = Lp + slab(i); 
end


%% now plot a single period of the Extended Basis

Nrepets = 1;

% nr of plot points
Npts = 10*round(Nrepets*Lp*10) + 1;
x = linspace(0,Nrepets*Lp,Npts)'; 
dx = x(2) - x(1) ;

V = hb*ones(Npts,1); 

offsets = zeros(length(slab),1);

for i =  2:length(slab)
        offsets(i) = offsets(i-1) + slab(i-1);
end

%Structure starts with a barrier! 
for i = 1:Npts
    period = ceil(x(i)/Lp);
    for j = 2:length(offsets)
        if((x(i)-(period-1)*Lp) > offsets(j-1) && (x(i)-(period-1)*Lp) <= offsets(j))
            if(mod(j-1,2) == 0) %if it is a well! 
                V(i) = 0; %        V(i) = V(i) -hb;
            end
            break
        end
    end
end

V = V+Bias*x;
%here we have a single period of the usual extended basis. 

%% read in the wavefunctions, the energies and the grid from the file.
xTB1 = xWF; TB_WF1 =wFunctions; TBEn = Energies; LpTB = xWF(end)-xWF(1); dxTB = xTB1(2)-xTB1(1); dLp = LpTB-Lp;

for i = 1:NrWF
    TB_WF1(:,i) = wFunctions(:,i)/sqrt(trapz(abs(wFunctions(:,i)).^2)*dxTB);
end

% now extend V,x and TB_WF2 to arbitrary number of  periods! 
Nper = 5; 
x = reshape(x,length(x),1);  V = reshape(V,length(V),1); 
x_fin = [x]; V_fin = [V]; 
for i = 1:Nper-1
    x_fin = [x_fin; (x+i*(Lp+dx))];
    V_fin = [V_fin; (V +i*(Lp+dx)*Bias)];
end

plot(x_fin,V_fin);
%%%%%%%%%%%%% extend V %%%%%%
extension_thicknesR = 2.4;  %nm 
dx = x_fin(2)-x_fin(1);
extension_pts = round((extension_thicknesR/dx)) + 1;

%pad the potential from the right! 
padxR = x_fin(end) +dx*(1:extension_pts);
% padVR = zeros(extension_pts,1);
padVR = V_fin(end) + [1:1:extension_pts].'*Bias*dx; 

% for i =2 : extension_pts
%     padVR(i) = padVR(i-1) +Bias*dx;
% end

V_fin = [V_fin;padVR];
x_fin = [x_fin ; padxR.'];
%%%%
central_per = 2;
Vtb = V_fin;
for i = 2:length(x_fin)
    period = ceil(x_fin(i)/Lp);
    if(period ~= central_per)
        Vtb(i) =  Vtb(i-1)+Bias*dx;
    end
end
figure;
% plot(x_fin,V_fin,'Linewidth',2.0,'color','k'); hold on; plot(x_fin,Vtb,'--r','Linewidth',2.0); xlim([x_fin(1),x_fin(end)])
hold on; 

%now interpolate the WFs onto the final grid! 
TB_WF2 = zeros(length(x_fin),NrWF,Nper); 
Ens  = zeros(NrWF,Nper); TB_WF2 = TB_WF2; 

B =[0,0,1;0,0.5,0;1,0,0;0,0.75,0.75;0.75,0,0.75;0.75,0.75,0;0.25,0.25,0.25;0,1,0;0.5,0.5,0.5;0.5,0.5,0;1,0.5,0];
A=B(1:NrWF,:);
set(gca,'ColorOrder',A);
for p =1:Nper
for i = 1:NrWF
        TB_WF2(:,i,p) = interp1(xTB1-dLp/2+(p-1)*(Lp),TB_WF1(:,i),x_fin,'linear',0);
        Ens(i,p) = Energies(i)+(p-1)*Lp*Bias;
%         plot(x_fin,abs((TB_WF2(:,i,p))).^2+Ens(i,p),'-','Linewidth',2);
end
end


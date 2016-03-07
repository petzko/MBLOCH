%%%% extract the IV characteristic of the laser

% width [um]
wum = 20; wcm = wum*1E-4;
%hegiht [um] 183 x 54.8 nm *1E-3; 
hum = 183*54.8*1E-3;
hcm = hum*1E-4; 
%length [um] 5mm x 1E3 
lum = 5*1E3; 
lcm = lum*1E-4; 

Biaslvls = [6.0 7.0 8.0 9.0 10.2 10.4 10.7 11.0 11.3 12.0 13.0 14.0]; %in kV/cm
%Biaslvls = [6.0 8.0 9.0 10.4 10.7 11.0 11.3 12.0 13.0 14.0]; %in kV/

nlvls = length(Biaslvls);
V = [] ;

for vidx = 1:nlvls
    V(vidx) = Biaslvls(vidx)*1E3*hcm;
end

maindir = pwd; 

%current density in A/m^2
A = wcm*1E-2*lcm*1E-2
J = []; 
for vidx = 1:nlvls  
    dir = [maindir '\' num2str(Biaslvls(vidx),'%.1f')];
    [err,time, Itmp] = currPETZ(dir,100,5E4,0); 
    J(vidx) = mean(Itmp);
end

I = J*A; 

plot(V,I);
xlabel('V (Volts)'); ylabel('I (Amperes)'); 





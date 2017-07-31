clear; clc; close all; 
TRT = 30; % ps
T1g = 10; % ps 
T2g = 0.5; % ps 
zULg  = 2; 
T2a = 0.5; %ps 

zULa = linspace(1,10,10); 
T1a = linspace(0.5,5,10); 
FML_AREAs = zeros(length(T1a),length(zULa)); 
HML_AREAs = zeros(length(T1a),length(zULa)); 
% 
for m=1:length(T1a)
    for n=1:length(zULa)
     [FML,HML] = new_condition(T1g,T1a(m),T2g,T2a,zULg,zULa(n),TRT);
     FML_AREAs(m,n) = FML; 
     HML_AREAs(m,n) = HML;
    end
end

save('optim_area_scan_T1g=10ps;T2g=0p5ps;T2a=0p5ps;zULg=2nm;TRT=30ps;zULa=varied;T1a=varied;'); 

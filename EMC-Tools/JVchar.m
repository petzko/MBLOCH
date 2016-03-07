

Biaslvls = [38 40 43 45 48 50 53]; %in kV/cm

nlvls = length(Biaslvls);
V = [] ;

for vidx = 1:nlvls
    V(vidx) = Biaslvls(vidx);
end

maindir = 'D:\Ezhovivan\Simulation\Screening\Mid-infra\New_order';  

%current density in A/m^2

J = []; 
for vidx = 1:nlvls  
    dir = [maindir '\' num2str(Biaslvls(vidx),'%d')];
    [err,time, Itmp] = currPETZ(dir,100,5E4,1); 
    J(vidx) = mean(Itmp);
end


plot(V,J,'-xr','Linewidth',2.0);
xlabel('V (kV/cm)'); ylabel('J(A/m^2)'); 





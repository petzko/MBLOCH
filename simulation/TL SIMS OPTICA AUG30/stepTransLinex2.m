function [ dat ] = stepTransLinex2(settings,dat)

%
shb11 = dat.r11p.*exp(2*1i*dat.x*dat.E0/dat.c);
shb33 = dat.r33p.*exp(2*1i*dat.x*dat.E0/dat.c);
shb22 = dat.r22p.*exp(2*1i*dat.x*dat.E0/dat.c);

shb13 = dat.r13p.*exp(2*1i*dat.x*dat.E0/dat.c) + dat.r13m.*exp(-2*1i*dat.x*dat.E0/dat.c);

%current density will be calculated with SHB included! 
r11 = dat.r110+2*real(shb11);     r33 = dat.r330+2*real(shb33);  r22 = dat.r220+2*real(shb22);
r13 = dat.r130+shb13;
    
dat.rates1 =   r11.*dat.W(:,dat.INJ,dat.DEPOP) +  r33.*dat.W(:,dat.ULL,dat.DEPOP)+...
    r22.*dat.W(:,dat.LLL,dat.DEPOP)+dat.rRES.*dat.W(:,dat.RES,dat.DEPOP);
dat.rates2 =   1i*dat.O13.*(r13-conj(r13));

dat.J_TL1 =  (Constants('q0')*(settings.Lp*1E-9)*dat.Ncarriers*1E12*dat.rates1)/1E6; %in A/mm^2
dat.J_TL2 =  -(Constants('q0')*(settings.Lp*1E-9)*dat.Ncarriers*1E12*dat.rates2)/1E6; %in A/mm^2

J_TL_old = dat.J_TL(1); 
dat.J_TL =dat.J_TL2;

%update current -> temporal grid 1,2,.. spatial grid 1/2,3/2 ...
% dat.i_TL(1:end) = dat.i_TL(1:end)*(1-dat.dt*dat.R_Au)-dat.Dcoeff*(dat.v_TL(2:end)-dat.v_TL(1:end-1));
dat.i_TL(1:end) = dat.i_TL(1:end)-dat.Dcoeff*(dat.v_TL(2:end)-dat.v_TL(1:end-1));

%update voltage -> temporal grid 1/2,3/2,5/2 ..., spatial grid 0,1,2 ...
dat.v_TL(2:end-1) = dat.v_TL(2:end-1)-dat.Ecoeff*(dat.i_TL(2:end)-dat.i_TL(1:end-1))-dat.Fcoeff*dat.J_TL(2:end-1);

%bdry conditions
%left side -> Kirchoff's law
v_1_old = dat.v_TL(1); 
dat.v_TL(1) = 1/(dat.width_mm*dat.dx/2/dat.Fcoeff+dat.height_mm/2/dat.Z_kV_A)*...
            (dat.Vs(dat.t)/dat.Z_kV_A-dat.i_TL(1)*dat.width_mm - ... %injection voltage and current 
            dat.height_mm*dat.width_mm*(dat.J_TL(1)+J_TL_old)/2 + ...% current density*Area! 
            (dat.width_mm*dat.dx/2/dat.Fcoeff-dat.height_mm/2/dat.Z_kV_A)*v_1_old);
        
%right side -> OC
dat.v_TL(end) = dat.v_TL(end)+2*dat.Fcoeff/dat.dx*dat.i_TL(end)-dat.Fcoeff*dat.J_TL(end);
% %

end


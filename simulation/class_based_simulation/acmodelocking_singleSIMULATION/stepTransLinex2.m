function [ dat ] = stepTransLinex2(params,dat,J_TL)

%update current -> temporal grid 1,2,.. spatial grid 1/2,3/2 ...
% dat.i_TL(1:end) = dat.i_TL(1:end)*(1-dat.dt*dat.R_Au)-dat.Dcoeff*(dat.v_TL(2:end)-dat.v_TL(1:end-1));
dat.i_TL(1:end) = dat.i_TL(1:end)-dat.Dcoeff*(dat.v_TL(2:end)-dat.v_TL(1:end-1));

%update voltage -> temporal grid 1/2,3/2,5/2 ..., spatial grid 0,1,2 ...
dat.v_TL(2:end-1) = dat.v_TL(2:end-1)-dat.Ecoeff*(dat.i_TL(2:end)-dat.i_TL(1:end-1))-dat.Fcoeff*J_TL(2:end-1);

%bdry conditions
%left side -> Kirchoff's law
% v_1_old = dat.v_TL(1); 
% dat.v_TL(1) = 1/(dat.width_mm*dat.dx/2/dat.Fcoeff+dat.height_mm/2/dat.Z_kV_A)*...
%             (dat.Vs(dat.t)/dat.Z_kV_A-dat.i_TL(1)*dat.width_mm - ... %injection voltage and current 
%             dat.height_mm*dat.width_mm*(dat.J_TL(1)+J_TL_old)/2 + ...% current density*Area! 
%             (dat.width_mm*dat.dx/2/dat.Fcoeff-dat.height_mm/2/dat.Z_kV_A)*v_1_old);
dat.v_TL(1) = dat.v0; 
%right side -> OC
dat.v_TL(end) = dat.v_TL(end)+2*dat.Fcoeff/params.dx*dat.i_TL(end)-dat.Fcoeff*J_TL(end);
% %

end


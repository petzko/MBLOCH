function Vx_dope = CalcVx_dope(Psi,E,dopx,rindexx,dx)

nx = length(Psi(:,1));
const = Constants();

charge_level = 4;

ncarriers = sum(dopx*dx); % total number of carriers per unit cross-sectional area
ncarriers

Psin = Psi(:,charge_level);

for ii = 1 : nx
    netchargedensity(ii) = const.qe * (ncarriers*Psin(ii)*Psin(ii)-dopx(ii));
end

efield_dope = [ ];

for ii = 1 : nx
    if ii == 1
        efield_dope(ii) = netchargedensity(ii)*dx;
    else
        efield_dope(ii) = efield_dope(ii-1) + netchargedensity(ii)*dx;
    end
end

Vx_dope = [ ];

for ii = 1 : nx
    if ii == 1
        Vx_dope(ii) = efield_dope(ii) * dx;
    else
        Vx_dope(ii) = Vx_dope(ii-1) + efield_dope(ii) * dx;
    end
end

for ii = 1 : nx
    Vx_dope(ii) = Vx_dope(ii)/(rindexx(ii)^2*const.eps0);
end

x = (0 : (nx-1))'*dx;
%Vx_test = Vx+Vx_dope;
figure(4)
plot(x,Vx_dope)
%figure(5)
%plot(x,efield_dope)
%figure(6)
%plot(x,dopx)
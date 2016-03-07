function [x_0,per] = centroid(Psi,x,Lp)
    %normalize 
    Psi = Psi/sqrt(trapz(x,abs(Psi).^2));
    x_0 = trapz(x,x.*abs(Psi).^2);
    per = floor(x_0/Lp)+1;
end
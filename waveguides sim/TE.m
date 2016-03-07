function erg=TE(bet)
global par; global M1; global k1; global k2;
beta=bet(1)+i*bet(2); k0=par(1); d=par(2); eps1=par(3); eps2=par(4); 
k1=sqrt(eps1*k0.^2-beta^2); k2=sqrt(eps2*k0.^2-beta^2); if(imag(k1)<0)k1=-k1; end; if(imag(k2)<0)k2=-k2; end; k3=k1;
%impl=(g2-g3)*(g2-g1)*exp(2*i*k2*d)-(g2+g3)*(g2+g1);
M1=[(1+k1/k2)/2 (1-k1/k2)/2;(1-k1/k2)/2 (1+k1/k2)/2]; M2=[exp(i*k2*d) 0;0 exp(-i*k2*d)]; M3=[(1+k2/k3)/2 (1-k2/k3)/2;(1-k2/k3)/2 (1+k2/k3)/2];
M=M3*M2*M1;
%erg=[real(impl) imag(impl)];
erg=[real(M(2,2)) imag(M(2,2))];
function erg=TM7(bet)
global k0; global b; global eps; global ig;
beta=bet(1)+i*bet(2);
k=sqrt(eps*k0.^2-beta^2);
k=k.*(1-2*(imag(k)<0));
g=k./eps;
for n=1:(length(eps)-1)
Mn(:,:,2*n-1)=[(1+g(n)/g(n+1))/2 (1-g(n)/g(n+1))/2;(1-g(n)/g(n+1))/2 (1+g(n)/g(n+1))/2];
if(n<(length(eps)-1)) Mn(:,:,2*n)=[exp(i*k(n+1)*b(n+1)) 0;0 exp(-i*k(n+1)*b(n+1))]; end;
end;
M=Mn(:,:,1);
for n=2:(2*length(eps)-3) M=Mn(:,:,n)*M; end;
%erg=[real(impl) imag(impl)];
erg=[real(M(2,2)) imag(M(2,2))];
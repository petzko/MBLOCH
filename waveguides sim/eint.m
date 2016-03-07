function erg=eint(M,k,d,sw);
%for sw=1 take the derivative
AB=M*[0;1]; kr=real(k); ki=imag(k);
if(sw==1) AB=i*k*AB; AB(2)=-AB(2); end;
erg=abs(AB(1))^2/2/ki*(1-exp(-2*d*ki))+abs(AB(2))^2/2/ki*(exp(2*d*ki)-1)+real(AB(1)*conj(AB(2))/i/kr*(exp(2*i*d*kr)-1));
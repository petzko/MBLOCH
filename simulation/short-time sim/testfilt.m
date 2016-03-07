x=randn(100,1)+1i*randn(100,1);
B=fir1(10,.1);

xfri=filter(B,1,real(x))+1i*filter(B,1,imag(x));
xf=filter(B,1,x); % Filter Real and Imag together

plot(real(xfri'));
hold on
plot(real(xf')) 
% plot([imag(xfri') imag(xf')]) % Plot Imag Parts


L = 1;
x = linspace(0,L,1000);
k = 10;
lambda = 2*pi/k;
A = 1; pi/2

D_2 = @(x) 19.5+ sin(k.*x+pi/2).^2;
I = @(x) A.^2*sin(k.*x).^2;
D_0 = @(x) 20 + 0.00001*sin(k.*x+pi/2).^2;

d0 = D_0(x); d2 = D_2(x);
subplot(2,1,1)
plot(x,d0,'r',x,d2,'g');
subplot(2,1,2);
plot(x,I(x),'b');

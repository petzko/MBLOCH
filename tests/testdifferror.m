
clear;
v = 1;
w = 16;
k = w/v;
a = -1;
b = 1;
% 
% f_x = @(x,t) exp(1i*(k*x-w*t)) + 1i*k*(x).*exp(1i*(k*x-w*t));
% f = @(x,t) (x + v*t).*exp(1i*(k*x-w*t));

f = @(x) sin(x);
f_x = @(x) cos(x);


orderend = 7; 
init_size = 64;
sizes = zeros(orderend,1);
sizes(1) = init_size;

t = 0.0;

for s = 2: orderend
    sizes(s) = sizes(s-1)*2;
end



err = 0;

dx = 0
for order = 1:orderend
    
    N = sizes(order);
%     x = -cos(pi*(0:(N-1))/(N-1))';
    [D_N,x] = cheb(N-1);
    dx(order) = x(2) - x(1);
    err(order) = max(abs(D_N*f(x) - f_x(x)));

end

display(err);
rho = 0;
for k = 2:length(err)
    
    rho(k-1) = log(err(k-1)/err(k))/log(dx(k-1)/dx(k));
    
end


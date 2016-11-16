close all;
c = Constants('c');
eps0 = Constants('eps0'); %F/m
q0 = Constants('q0'); % C
hbar = Constants('hbar');  % Js
deg = 4.6e-9; %m

Olap = .95;   
nthz = 3.6;
Ose =  0.243e-3*q0/hbar;
Ose = 0.0587*2*pi*1e12;
Ose = 0.0563*2*pi*1e12;

epsilon = -0*1e-3*q0/hbar;
% wc = 20.85e-3*q0/hbar; %1/s; 
wc =  5.0433*1e12*2*pi
%slow light
N = 7.8946e21;%in m^-3; 
r= -1

Op = N*Olap*(q0*deg)^2/(2*nthz^2*eps0*hbar);
ng_est = 1+abs(r)*Op*wc/Ose^2

%%%% two level gain coeff %%%% 
T2 = 2.35e-12;
Ngain = 2.1e17*(100)^3;
zUL_GAIN = 4.5e-9;
GAIN = T2*wc*Ngain*(q0*zUL_GAIN)^2/(hbar*eps0*c*nthz)
GAIN*0.01;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

df = 1e7; 
f = [-100000:100000]*df; 

w=f*2*pi;
dw = w(2)-w(1);

dL = .1e-3; %.1 mm
c = Constants('c');
xl = [-.4e12,.4e12];
colors = {'b','r'};
for j =0:1 %switch dephasing on/off
    col = colors{j+1};
    Gsg = j*0.01801e12; Geg = 0.2038e12 ;
%     Gsg = j*0.0357e12; Geg = 0.2163e12 ;
%     Gsg = j*1e12; Geg = 1e12 ;
    n_ = @(w_) nthz*(1+r*Op*(w_-epsilon +1i*Gsg)./((w_-epsilon +1i*Gsg).*(w_ +1i*Geg)-Ose^2));
    n_im =@(w_) -r*nthz*Op*(Geg*(Gsg^2+w_.^2)+Ose^2*Gsg)./((w_.^2-Geg*Gsg-Ose^2).^2+w_.^2*(Geg+Gsg).^2);
    n_g = @(w_) real(n_(w_(1:end-1)))+wc*real(diff(n_(w_))/dw);
    ng = nthz*(1-r*wc*Op*(Ose^2-Gsg^2)/(Ose^2+Gsg*Geg)^2);
    T_w = @(w_) exp(-2.*(w_+wc)/c.*imag(n_(w_))*dL);

    subplot(4,1,1); hold on;
    plot(f,real(n_(w)),col); ylabel('n_{re}');xlabel('\omega/2\pi')
    xlim(xl);
    
    subplot(4,1,2); hold on;
    alpha = (2*imag(n_(w)).*(w+wc)/c)*0.01; %per cm 
    [a_min,frq] = min(alpha( f> -1000*df & f<1000*df));
    display(['Min absorption = ' num2str(a_min) ' at freq ' num2str(frq)]);  

    plot(f,alpha); ylabel('\alpha(\omega) 1/cm')
    xlim(xl);
    
    subplot(4,1,3);hold on;
    ng_w = n_g(w);
    [ng_max,frq] = max(ng_w( f> -1000*df & f<1000*df)); 
    display(['Max delay = ' num2str(ng_max) ' at freq ' num2str(frq)]);  
    plot(f(1:end-1),ng_w,col);
    ylabel('n_{g}');xlabel('\omega/2\pi')
    xlim(xl);
    
    subplot(4,1,4);hold on;
    plot(f,T_w(w),col);ylabel('Transmittance');xlabel('\omega/2\pi')
    xlim(xl);
    
end

display(['L/rOp: ' num2str(((nthz-1)/nthz)*Ose^2/wc/(r*Op) ) ' R/rOp: ' num2str(Ose^2/wc/(r*Op)) '; rOp = ' num2str(r*Op)]);

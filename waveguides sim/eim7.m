function erg=eim7(f,chi,T,sw,c1)
%a=25e-6; ymin=2.0950e-6; ymax=3.8860e-6; Seff=trapz(e1(1,:),abs(e1(2,:)).^2)*trapz(e1(1,:),abs(e2(2,:)).^2)*trapz(e1(1,:),abs(e3(2,:)).^2)/abs(trapz(e1(1,:),e1(2,:).*e2(2,:).*e3(2,:).*(e1(1,:)>ymin).*(e1(1,:)<ymax)))^2*a;
eps0=8.8541878176e-12; e0=1.60217646e-19; hbar=1.05457168e-34; c0=3e8; me=9.10938188e-31;
global k0; global b; global eps; global ig; global aw; global Gam1; global T12;
cx=0.52; %Al_{cx}Ga_{1-cx}As
mo=0; %mo is mode index in x direction
if(~exist('T')) T=300; end;
if(~exist('sw')) sw=1; end;

t_m=0.06e-12; t_c=0.1e-12; t_w=0.5e-12; %koh05: t_m for metal; t_c for contact layer (e.g., N=5e24); t_w for weakly doped semiconductor (e.g., N=5.5e21)
w=2*pi*f; k0=w/c0; wcm=f/3e10; %wcm is frequency in cm^-1
%GaAs=1; AlAs=2; InAs=3; InP=4; In_{0.53}Ga_{0.47}As=5; Al_{cx}Ga_{1-cx}As=6; gain medium=7; Au=8; Cu=9; Ti=10;
epsm=[NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];
if(f>2.7e13) %used interpolation formulas good for lambda<11µm
%IR
la=c0/f*1e6; 
%GaAs
C0=3.5; C1=7.4969; C2=1.9347; l1=0.4082; l2=37.17;
epsm(1)=C0+C1*la^2/(la^2-l1^2)+C2*la^2/(la^2-l2^2);
%InAs
C0=11.1; C1=0.71; C2=2.75; l1=2.551; l2=45.66;
epsm(3)=C0+C1*la^2/(la^2-l1^2)+C2*la^2/(la^2-l2^2);
%InP
C0=7.255; C1=2.316; C2=2.765; l1=0.6263; l2=32.935;
epsm(4)=C0+C1*la^2/(la^2-l1^2)+C2*la^2/(la^2-l2^2);
%InGaAs
% C0=10.65; C1=1.35; C2=4258; l1=1.24; l2=1020;
% epsm(5)=C0+C1*la^2/(la^2-l1^2)+C2*la^2/(la^2-l2^2);
A0=7.4326; B0=-0.0077897; C0=9.6836;
E=hbar*2*pi*3e8./(la*1e-6)/e0; E0=0.737; D0=0.354; x0=E/E0; xos=E/(E0+D0);
F0=(2-sqrt(1+x0)-sqrt(1-x0).*(x0<=1))./x0.^2; Fox=(2-sqrt(1+xos)-sqrt(1-xos).*(xos<=1))./xos.^2;
epsm(5)=A0*(F0+0.5*(E0./(E0+D0)).^1.5*Fox)+B0./E.^2+C0;
%Al_{cx}Ga_{1-cx}As
C0=1; C1=(36.1-2.45*cx)/(3.65+0.871*cx+0.179*cx^2); C2=0; l1=1.241/(3.65+0.871*cx+0.179*cx^2); l2=0;
epsm(6)=C0+C1*la^2/(la^2-l1^2)+C2*la^2/(la^2-l2^2);
%Au
%N_Au=5.9e28; wp_Au=sqrt(e0^2.*N_Au./eps0./me); epsm(8)=1-wp.^2./(w.^2+i*w./t_m); %Drude model for metals
wt=216; wpf=72500; %parameters from AO 22, 1099
epsm(8)=-wpf^2/(wcm^2+i*wcm*wt);
%Cu
wt=278; wpf=63800; %parameters from AO 22, 1099
epsm(9)=-wpf^2/(wcm^2+i*wcm*wt);
%Ti
wt=732.7; wpf=22178; epsi=i*51.62; %parameters fitted from AO 22, 1099
epsm(10)=-wpf^2/(wcm^2+i*wcm*wt);
end;

if(f<10e12)
%THz
%InP
ELO=43e-3; ETO=38.1e-3; G=0.166e-3; epsi=9.61; %InP
E=(2*pi*hbar*f)/e0; ef=epsi*(1+(ELO^2-ETO^2)/(ETO^2-E^2-i*G*E));
epsm(4)=ef;
%In_0.53Ga_0.47As=In_1-xGa_xAs
epsix=11.6; eps0x=13.9; eps0x0=15.15; eps0x1=12.95;
ETOx0=28.3e-3; ETOx1=31.6e-3; G=0.0395e-3;
E=(2*pi*hbar*f)/e0; ef=epsix+ETOx0^2*(eps0x-eps0x1)/(ETOx0^2-E^2-i*G*E)+ETOx1^2*(eps0x1-epsix)/(ETOx1^2-E^2-i*G*E);
epsm(5)=ef;
%Al_{cx}Ga_{1-cx}As
epsix=10.89-2.73*cx; eps0x=13.18-3.12*cx; eps0x0=13.18; eps0x1=10.06;
ETOx0=(33.29-0.64*cx-1.16*cx^2)*1e-3; ETOx1=(44.63+0.55*cx-0.30*cx^2)*1e-3; G=0.2e-3; %values of 0.2 for AlAs and 0.23 or 0.3 for GaAs are found; however, the result is insensitive to G
E=(2*pi*hbar*f)/e0; ef=epsix+ETOx0^2*(eps0x-eps0x1)/(ETOx0^2-E^2-i*G*E)+ETOx1^2*(eps0x1-epsix)/(ETOx1^2-E^2-i*G*E);
epsm(6)=ef;
%Au
epsm(8)=(0.1516*(f/1e12).^2-0.8827*(f/1e12)-8.6609)*1e4+i*(0.1072*(f/1e12).^2-1.5998*(f/1e12)+5.8608)*1e5; %1..10THz, fitted to data in AO 22, 1099
%Ti
wt=233.2; wpf=20780; epsi=i*120.27; %parameters fitted from Palik, Handbook of Optical Constants of Solids, Band 3, p. 248-249
epsm(10)=-wpf^2/(wcm^2+i*wcm*wt);
end;

%Mobility
%GaAs=1; AlAs=2; InAs=3; InP=4; In_{0.53}Ga_{0.47}As=5; Al_{cx}Ga_{1-cx}As=6;
mumin=[500 10 1000 400 300 NaN NaN NaN NaN NaN]*1e-4; mumax=[9400 400 34000 5200 14000 NaN NaN NaN NaN NaN]*1e-4; Neff=[6.0e16 5.46e17 1.1e18 3.0e17 1.3e17 NaN NaN NaN NaN NaN]*1e6;
t1=[2.1 2.1 1.57 2.0 1.59 NaN NaN NaN NaN NaN]; t2=[3.0 3.0 3.0 3.25 3.68 NaN NaN NaN NaN NaN]; lam=[0.394 1 0.32 0.47 0.48 NaN NaN NaN NaN NaN];

%Effective mass
mm=[NaN NaN NaN 0.077 0.039 NaN NaN NaN NaN NaN];

%Parameters Razeghi APL 99, 131106
eps1=epsm(8); global dw; a=dw; %a=60e-6; %waveguide width
vm=[4 4 7 5 7 4 4 8]; %layer sequence
%Gain medium effective mass
cw=0.70; epsm(7)=1/(cw/epsm(5)+(1-cw)/epsm(6)); %cw=d_well/d_period is approximately 0.70 for both gain media
mm(7)=0.048; %Interpolation between InGaAs and InAlAs contributions yields this result for both the BTC and RP structure
mumin(7)=mumin(5); mumax(7)=mumax(5); t1(7)=t1(5); t2(7)=t2(5); lam(7)=lam(5); Neff(7)=Neff(5);
%end gain medium effective massN=[1.5e23 1.5e22 5e22 1e22 5e22 1.5e22 5e24 0];
meff=mm(vm);
N=[1.5e23 1.5e22 0 1e22 0 1.5e22 5e24 0]; %Metal and doping in gain medium (N=5e22) is treated separately
ig=[3 5]; %regions which are gain media
b=[NaN 5e-6 30*597e-10 100e-9 30*665e-10 3.5e-6 200e-9 NaN];
%g0=[0 0 -30000 0 -30000 0 0 0];
if(~exist('c1')) c1=1000; end;
b(1)=c1*b(ig(1)); b(end)=c1*b(ig(1)); %for plotting
%N=fliplr(N); tau=fliplr(tau); eps_u=fliplr(eps_u); meff=fliplr(meff); b=fliplr(b);
%end parameters

% %Parameters Belkin 2008 APL 92, 201101
% a=25e-6; %waveguide width
% vm=[4 7 5 7 5 4 4 10 8]; %layer sequence
% %Gain medium effective mass
% cw=0.70; epsm(7)=1/(cw/epsm(5)+(1-cw)/epsm(6)); %cw=d_well/d_period is approximately 0.70 for both gain media
% mm(7)=0.048; %Interpolation between InGaAs and InAlAs contributions yields this result for both the BTC and RP structure
% mumin(7)=mumin(5); mumax(7)=mumax(5); t1(7)=t1(5); t2(7)=t2(5); lam(7)=lam(5); Neff(7)=Neff(5);
% %end gain medium effective mass
% meff=mm(vm);
% N=[9e22 0 3e22 0 3e22 5e22 5e24 0 0]; %Metal and doping in gain medium (N=5e22) is treated separately
% ig=4; %ig=[2 4]; %regions which are gain media; used e.g. for calculating overlap factor
% b=[NaN 30*665e-10 100e-9 30*597e-10 50e-9 3.5e-6 200e-9 20e-9 NaN];
% b1=4*b(ig(1)); bend=0.05*b(ig(1)); %for plotting
% %end parameters

% %Parameters OL 32, 2840 (2007) QCL, but with a=100e-6 instead of a=80e-6 and f=4.5e12 instead of f=4.1e12;
% mo=0; n0=3.8; a=100e-6; %mo is mode index in x direction
% N=[5.9e28 5e24 5.5e21 5e24 5.9e28];
% tau=[0.06e-12 0.1e-12 0.5e-12 0.1e-12 0.06e-12];
% eps_u=[1 12.96 12.96 12.96 1]; %undoped eps
% ig=[3]; %regions which are gain media
% meff=[1 0.067 0.067 0.067 1];
% b=[NaN 100e-9 10e-6 50e-9 NaN];
% %end parameters

% %Parameters Wil EL QCL, but with a=100e-6 instead of a=98e-6 and f=4.5e12
% mo=0; n0=3.8; a=100e-6; %mo is mode index in x direction
% N=[5.9e28 0 5.5e21 3e24 0e24];
% tau=[0.06e-12 0.5e-12 0.5e-12 0.1e-12 0.5e-12];
% eps_u=[1 12.96 12.96 12.96 12.96];
% meff=[1 0.067 0.067 0.067 0.067];
% b=[NaN 100e-9 10e-6 400e-9 NaN];
% ig=[3]; %regions which are gain media
% %end parameters

%tau=[t_w t_w t_w t_w t_w t_w t_c t_m];
%tau=[t_c t_c t_c t_c t_c t_c t_c t_c]; %see Nat. Photonics supplement
%muL=0.57; N0=1e23; tau=muL./(1+sqrt(N/N0)).*meff*me/e0; %Electron. Lett. 10, 259 (for InP at room temperature)
mu=mumin(vm)+(mumax(vm).*(300/T).^t1(vm)-mumin(vm))./(1+(N./Neff(vm).*(300/T).^t2(vm)).^lam(vm)); tau=mu.*meff*me/e0;
eps_u=epsm(vm); %undoped eps

if(~exist('chi')) chi=0*eps_u; end;
wp=sqrt(e0^2.*N./eps0./me./meff);
eps_drude=-wp.^2./(w.^2+i*w./tau);
for n=1:length(eps_drude) if(N(n)==0) eps_drude(n)=0; end; end;
eps=eps_u+eps_drude+chi; %chi=-i*ng.*g0/k0

%myfun = @(x)modeeq(x,k0,b,eps1,eps2)
options=optimset('TolFun',1e-20);
%1D problem
[bet,fval] = fsolve(@TM7,[0.99*real(sqrt(eps_u(ig(1))))*k0 0],options); beta=bet(1)+i*bet(2);
aw=2*imag(beta)/100
neff=beta/k0
%Overlap
k=sqrt(eps*k0.^2-beta^2);
k=k.*sign(imag(k));
g=k./eps;
ki=imag(k(1));
K=length(eps);
for n=1:(K-1)
    Mn(:,:,2*n-1)=[(1+g(n)/g(n+1))/2 (1-g(n)/g(n+1))/2;(1-g(n)/g(n+1))/2 (1+g(n)/g(n+1))/2];
    if(n<(K-1)) Mn(:,:,2*n)=[exp(i*k(n+1)*b(n+1)) 0;0 exp(-i*k(n+1)*b(n+1))]; end;
    %Overlap factor
    if(n==1)
        Ey2=abs(beta/eps(1))^2/2/ki; Ez2=abs(k(1)/eps(1))^2/2/ki; E2=Ey2+Ez2; E2g=zeros(1,length(ig));
    else
        if(n==2) Mh=Mn(:,:,1); else Mh=Mn(:,:,2*n-3)*Mn(:,:,2*n-4)*Mh; end;
        Ey2=abs(beta/eps(n))^2*eint(Mh,k(n),b(n),0); Ez2=abs(1/eps(n))^2*eint(Mh,k(n),b(n),1); E2=E2+Ey2+Ez2;
        if(sum(ig==n)) [dummy,ind]=max(ig==n); E2g(ind)=E2g(ind)+Ey2; end;
    end;
    %%%%%%%%%%%%%%%%%
end;
ind=[];
ki=imag(k(K)); Mh=Mn(:,:,2*K-3)*Mn(:,:,2*K-4)*Mh*[0;1]; Ey2=abs(Mh(1)*beta/eps(K))^2/2/ki; Ez2=abs(Mh(1)*k(K)/eps(K))^2/2/ki; E2=E2+Ey2+Ez2;
Gam1=E2g/E2

% v=[1e-9 (0.01:0.01:1)];
% y(1,:)=-b1*(1:(-0.01):0); Hy(1,:)=exp(-i*k(1)*y(1,:));
% y(2,:)=v*b(2); Hh=Mn(:,:,1)*[0;1]; Hy(2,:)=Hh(1)*exp(i*k(2)*y(2,:))+Hh(2)*exp(-i*k(2)*y(2,:));
% for n=3:(K-1)
%     y(n,:)=v*b(n); Hh=Mn(:,:,2*n-3)*Mn(:,:,2*n-4)*Hh; Hy(n,:)=Hh(1)*exp(i*k(n)*y(n,:))+Hh(2)*exp(-i*k(n)*y(n,:));
% end;
% n=K; y(n,:)=bend*v; Hh=Mn(:,:,2*n-3)*Mn(:,:,2*n-4)*Hh; Hy(n,:)=Hh(1)*exp(i*k(n)*y(n,:))+Hh(2)*exp(-i*k(n)*y(n,:));
% for n=2:K y(n,:)=y(n,:)+y(n-1,end); end;
% yv=reshape(y',1,[]); Hyv=reshape(Hy',1,[]); Hyv=Hyv/sqrt(trapz(abs(Hyv).^2,yv));
% figure; plot(yv,abs(Hyv).^2);
% hold on;

yif(1)=0;
for n=2:K yif(n)=yif(n-1)+b(n); end;
ymin=-b(1); ymax=yif(end); Nf=2^16; y=ymin+(0:(Nf-1))/(Nf-1)*(ymax-ymin);
yh=y; ind=find(yh<=0); Hy(ind)=exp(-i*k(1)*yh(ind)); epsy(ind)=eps(1);
yh=y; ind=find((yh>0).*(yh<=b(2))); Hh=Mn(:,:,1)*[0;1]; Hy(ind)=(Hh(1)*exp(i*k(2)*yh(ind))+Hh(2)*exp(-i*k(2)*yh(ind))); epsy(ind)=eps(2);
for n=3:(K-1)
    yh=y-yif(n-1); ind=find((yh>0).*(yh<=b(n))); Hh=Mn(:,:,2*n-3)*Mn(:,:,2*n-4)*Hh; Hy(ind)=(Hh(1)*exp(i*k(n)*yh(ind))+Hh(2)*exp(-i*k(n)*yh(ind))); epsy(ind)=eps(n);
end;
n=K; yh=y-yif(n-1); ind=find(yh>0); Hh=Mn(:,:,2*n-3)*Mn(:,:,2*n-4)*Hh; Hy(ind)=Hh(1)*exp(i*k(n)*yh(ind)); epsy(ind)=eps(n);

% %test
% epsy=0*y+1; beta=1*k0;
% si=1e-4; Hy=exp(-y.^2/2/si^2);
% %%%%%%%%%%

if(sw>0)
    figure; subplot(3,1,1); plot(y,abs(Hy).^2); hold on;
for n=ig
    plot([yif(n-1) yif(n)],[0 0],'linewidth',3);
end;
end;
Ex=-beta/eps0./epsy/w.*Hy; erg=[y;Ex];

%Transmission coefficient from IEEE JQE 10, 809 with corrections described in resonator.tex; change Nf and y interval to assure convergence
ky=(-Nf/2:(Nf/2-1))*2*pi/(y(end)-y(1));

norm=abs(beta)/k0*trapz(y,(real(beta)*real(epsy)+imag(beta)*imag(epsy))./abs(epsy).^2.*abs(Hy).^2);
Hy=Hy/sqrt(norm);
Phi=(y(2)-y(1))/2/pi*fftshift(fft(fftshift(Hy))); Phis=(y(2)-y(1))/2/pi*fftshift(fft(fftshift(Hy./epsy)));
Hfh=sqrt(1-ky.^2/k0^2).*Phi.*Phis./(k0*sqrt(1-ky.^2/k0^2).*Phi+beta*Phis); ind=find((ky>=-k0).*(ky<=k0)); Hf=0*ky; Hf(ind)=Hfh(ind);
theta=asin(ky/k0);
if(sw>0)
subplot(3,1,2); plot(ky,abs(Phi).^2);
subplot(3,1,3); plot(theta(ind),abs(Hf(ind)).^2);
end;
T12=8*pi*k0*abs(beta)^3*trapz(theta(ind),abs(Hf(ind)).^2)
%T12=8*pi*abs(beta)^3*trapz(ky,abs(Hf).^2./sqrt(1-ky.^2/k0^2))
%
if(abs(sw)==2)
%EIM
global M1; global k1; global k2; global par; global Eyx; %for routine Seff; only field in gain medium, not in the side walls
v=[1e-9 (0.01:0.01:1)];
par=[k0,a,eps1,real(neff^2)];
h=0;
for m=1:10000
    h_old=h; hv=TE(m/1000*[k0 0]); h=hv(1);
    if(h*h_old<0)
    x2=v*a; Hx2=M1(1,2)*exp(i*k2*x2)+M1(2,2)*exp(-i*k2*x2); if(sum(Hx2(1:(end-1)).*Hx2(2:end)<0)==mo)break; end;
    end;
end;
if(m==10000)disp('required mode does not exist'); end;
par=[k0,a,eps1,neff^2];
[bet,fval] = fsolve(@TE,m/1000*[k0 0],options); beta2=bet(1)+i*bet(2);
aw=2*imag(beta2)/100
neff2=beta2/k0
%Overlap
k0=par(1); d=par(2); eps1=par(3); eps2=par(4); eps3=eps1;
k1=sqrt(eps1*k0.^2-beta2^2); k2=sqrt(eps2*k0.^2-beta2^2); k1=k1*sign(imag(k1)); k2=k2*sign(imag(k2)); k3=k1;
%impl=(g2-g3)*(g2-g1)*exp(2*i*k2*d)-(g2+g3)*(g2+g1);
M1=[(1+k1/k2)/2 (1-k1/k2)/2;(1-k1/k2)/2 (1+k1/k2)/2]; M2=[exp(i*k2*d) 0;0 exp(-i*k2*d)]; M3=[(1+k2/k3)/2 (1-k2/k3)/2;(1-k2/k3)/2 (1+k2/k3)/2];

x1=-0.3*(1:(-0.01):0)*a; Hx1=exp(-i*k1*x1);
x2=v*a; Hh=M1*[0;1]; Hx2=Hh(1)*exp(i*k2*x2)+Hh(2)*exp(-i*k2*x2);
x3=0.3*v*a; Hh=M3*M2*Hh; Hx3=Hh(1)*exp(i*k3*x3);
x=[x1 x2 x2(end)+x3]; Hx=[Hx1 Hx2 Hx3]; Ey=-beta2/eps0/w*[Hx1/eps1 Hx2/neff^2 Hx3/eps1];
if(sw>0) figure; plot(x,abs(Hx).^2/max(abs(Hx).^2),'b'); end;
% [X,Y]=meshgrid(x,y);
% [HX,HY]=meshgrid(Hx,Hy);
% %mesh(X,Y,abs(HX.*HY).^2)/max((max(abs(HX.*HY).^2)));
% xi=min(x)+[0:0.002:1]*(max(x)-min(x)); yi=min(y)+[0:0.02:1]*(max(y)-min(y));
% for n=2:length(x) if x(n)==x(n-1) x(n)=x(n-1)+1e-20; end; end;
% for n=2:length(y) if y(n)==y(n-1) y(n)=y(n-1)+1e-20; end; end;
% Hxi=interp1(x,Hx,xi);
% Hyi=interp1(y,Hy,yi); [HXi,HYi]=meshgrid(Hxi,Hyi);
% image(xi*1e6,yi*1e6,abs(HXi.*HYi).^2)/max((max(abs(HXi.*HYi).^2)));

ki=imag(k1); Ey2=1/2/ki; E2=Ey2;
Ey2=eint(M1,k2,a,0); E2=E2+Ey2; E2g=Ey2;
ki=imag(k3); AB=M3*M2*M1*[0;1]; Ey2=abs(AB(1))^2/2/ki; E2=E2+Ey2;
Gam2=E2g/E2
Gam=Gam1*Gam2
Eyx=[x2;-beta2/eps0/neff^2/w*Hx2/sqrt(trapz(x,abs(Ey).^2))];
%Transmissivity
global par2; global kyHf;
kyHf=[ky(ind);Phi(ind);Phis(ind)]; par2=[k0,a,beta2];
disp([length(ind),theta(ind(end))]);
tol=1e-12;
%T12=16*pi^2*k0^3*abs(beta2)^3*dblquad('kernel',0,theta(ind(end)),-pi,pi,tol)
%T12=16*pi^2*abs(beta2)^3*dblquad('kernel2',-k0,k0,-k0,k0,tol)
kx=(-1:0.01:1)*k0;
for n=1:length(kx) T12h(n)=trapz(ky(ind),kernel2(kx(n),ky(ind))); end;
T12=16*pi^2*abs(beta2)^3*trapz(kx,T12h)

% aw+Gam*g/100
% Kp=abs(trapz(x,abs(Hx).^2)*trapz(y,abs(Hy).^2)/trapz(x,Hx.^2)/trapz(y,Hy.^2))^2
end;
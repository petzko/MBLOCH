clear; clc; close all;
T1g = 15;
T1a = 15;

T2a = .5;
T2g = .5;

zULa = 8;
zULg = 2;

TRT = 30;


E0 = 3.5*2*pi;
nTHz = 3.6;
Ncarriers_cm_g = 1.9e15;
Ncarriers_cm_a = 3.7e14;
linear_loss = 0
Lg = 2.0;
La = .5;




mu_g = zULg*1E-9*Constants('q0');
mu_a = zULa*1E-9*Constants('q0');

sigma_g = T2g*(1e-12)*E0*1E12*mu_g^2/ ...
    (Constants('eps0')*nTHz*Constants('c')*Constants('hbar'));
sigma_a = T2a*(1e-12)*E0*1E12*mu_a^2/ ...
    (Constants('eps0')*nTHz*Constants('c')*Constants('hbar'));

Ncarriers_a = Ncarriers_cm_a*1e6;
Ncarriers_g = Ncarriers_cm_g*1e6;

Ga = sigma_a*Ncarriers_a*La*1e-3;
Gg = sigma_g*Ncarriers_g*Lg*1e-3;

powerloss = 2*linear_loss*100;
d_th =  powerloss/(sigma_g*Ncarriers_g)+Ga/Gg;
p = 1.1;

gth = d_th*sigma_g*Ncarriers_g*Lg*1e-3;
g0_ = p*gth
q0_ = -1*sigma_a*Ncarriers_a*La*1e-3;


%%
g0 = [p*gth:0.1:2.0*gth];
q0 = ([q0_:-.1:2*q0_]);
% q0 = q0_; 

Gamma = T2a/T2g;
s = (zULa/zULg).^2;
s = Gamma*s;

R =1.
kappa = R*exp(-powerloss*Lg*1e-3);

res_1 = zeros(length(g0),length(q0));
res_2 = zeros(length(g0),length(q0));
DPs = zeros(length(g0),length(q0));
G1s = zeros(length(g0),length(q0));
G2s = zeros(length(g0),length(q0));
Q1s = zeros(length(g0),length(q0));
Q2s = zeros(length(g0),length(q0));

options = optimoptions('fsolve','Display','none');
options.MaxIter = 1000;
options.MaxFunEvals = 1000;
options.Algorithm = 'levenberg-marquardt';
% options.TolFun = @(F) sum(abs(F));

tol = 1e-7;
MAXREP = 500; 

for m = 1:length(g0)
    m
    for n=1:length(q0)
%         n
        X0 = [rand(1,2), -rand(1,2), rand(1)];
        [X,fval,exitflag,output] = fsolve('background_stability',X0,options,s,TRT,T1g,T1a,g0(m),q0(n),kappa);
        ctr = 1;
        while ((sum(abs(imag(X))) > tol) || (X(end))<tol || exitflag <= 0) && (ctr < MAXREP)
            X0 = [rand(1,2), rand(1,2), 100*rand(1)];
            [X,fval,exitflag,output] = fsolve('background_stability',X0,options,s,TRT,T1g,T1a,g0(m),q0(n),kappa);
            if mod(ctr,50) == 0
                display(['repsol: ' num2str(ctr) ]);
            end
            ctr = ctr + 1;
        end
        X = real(X);
        
        G1 = X(1); G2 = X(2); Q1 = X(3); Q2 = X(4); DP =X(5);
        DP
        DPs(m,n) = DP;
        G1s(m,n) = G1; 
        G2s(m,n) = G2; 
        Q1s(m,n) = Q1; 
        Q2s(m,n) = Q2; 
        
        
        res_1(m,n) = G1+Q1+log(kappa) < 0;
        res_2(m,n) = G2+Q2+log(kappa) < 0;
    end
end
res = (res_1+0).*(res_2+0);
[Gs,Qs] = meshgrid(g0,q0);
subplot(2,1,1);
surf(Gs,Qs,res');
subplot(2,1,2);

DPs = DPs.*((DPs>0)+0.);
surf(Gs,Qs,DPs');


%%

% mu = 1e-2;
% T = mu/10;
%
% DPvals = linspace(-1.,1.0,10000);
%
% fermi_ = 1./(exp((DPvals-mu)/T)+1);
% figure;
% plot(DPvals/mu, fermi_);


% DP0 = 10;
% options = optimoptions('fsolve','Display','none','PlotFcn',@optimplotx );
% [DP,val] = fsolve('pulse_area_steady_state',DP0,options,s,kappa);

% DPvals = linspace(-10,10,1000);
% plot(DPvals,DPvals);
% hold on
% s = 5;
% f_ = kappa/s;
% % plot(DPvals,f_*(log(1+(1+exp(DPvals)).^s)- log(1+2^s)));



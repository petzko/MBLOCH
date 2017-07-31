clear; clc; close all;
T1g = 10;
T1a = 3.5;

gamma = T1a/T1g;
T2g = .5;
T2a = .5;


zULa = 10;
zULg = 2;

s = zULa^2*T1a*T2a/(zULg^2*T1g*T2g);


TRT = 30;
T = TRT/T1a;

E0 = 3.5*2*pi;
nTHz = 3.6;
Ncarriers_cm_g = 1.9e15;
Ncarriers_cm_a = 3.7e14;
linear_loss = 5;
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
g0 = gamma*linspace(gth,3*gth,30);
q0 = linspace(q0_,3*q0_,30); 

R = 1.
kappa = R*exp(-powerloss*Lg*1e-3);

% solutions for fundamental mode locking
DPs_T1 = zeros(length(g0),length(q0));
G1s_T1 = zeros(length(g0),length(q0));
G2s_T1 = zeros(length(g0),length(q0));
Q1s_T1 = zeros(length(g0),length(q0));
Q2s_T1 = zeros(length(g0),length(q0));


% solution for second harmonic mode locking
DPs_T2 = zeros(length(g0),length(q0));
G1s_T2 = zeros(length(g0),length(q0));
G2s_T2 = zeros(length(g0),length(q0));
Q1s_T2 = zeros(length(g0),length(q0));
Q2s_T2 = zeros(length(g0),length(q0));



options = optimoptions('fsolve','Display','none');
options.MaxIter = 1000;
options.TolX = 1e-12;
options.TolFun = 1e-12;
options.MaxFunEvals = 1000;
options.Algorithm = 'levenberg-marquardt';
% options.TolFun = @(F) sum(abs(F));

tol = 1e-10;
MAXREP = 500; 

for m = 1:length(g0)
    m
    for n=1:length(q0)
        
        % fundamental mode locking
        X_T1 = get_reasonable_solution(options,s,T,g0(m),q0(n),kappa,gamma,tol,MAXREP); 
        G1 = X_T1(1); G2 = X_T1(2); Q1 = X_T1(3); Q2 = X_T1(4); DP =X_T1(5);
        DPs_T1(m,n) = DP;   G1s_T1(m,n) = G1;    G2s_T1(m,n) = G2; 
        Q1s_T1(m,n) = Q1;   Q2s_T1(m,n) = Q2; 
        % second harmonic mode locking 
        X_T2 = get_reasonable_solution(options,s,T/2,g0(m),q0(n),kappa,gamma,tol,MAXREP);  
        G1 = X_T2(1); G2 = X_T2(2); Q1 = X_T2(3); Q2 = X_T2(4); DP =X_T2(5);
        DPs_T2(m,n) = DP;   G1s_T2(m,n) = G1;    G2s_T2(m,n) = G2; 
        Q1s_T2(m,n) = Q1;   Q2s_T2(m,n) = Q2; 
       
    end
end
%%
[Gs,Qs] = meshgrid(g0/gamma,q0);

FUNDAMENTAL_REGION = (G1s_T1'+Q1s_T1' + log(kappa) <= 0 & G2s_T1'+Q2s_T1' + log(kappa) <= 0) +0. ; 
SECOND_HARMONIC_REGION = (G1s_T2'+Q1s_T2' + log(kappa) <= 0 & G2s_T2'+Q2s_T2' + log(kappa) <= 0) +0.; 
% close all; 
THRESHOLD_REGION =  Gs + Qs + log(kappa);
THRESHOLD_REGION = THRESHOLD_REGION > 0;  
dTHRESHOLD = (diff(THRESHOLD_REGION) == -1) +0; 

dFUNDAMENTAL = (diff(FUNDAMENTAL_REGION,1) == 1)+0.;
dSECOND_HARMONIC = (diff(SECOND_HARMONIC_REGION,1) == 1) + 0; 


figure; 
contour(Gs(2:end,:),Qs(2:end,:),dTHRESHOLD,[1],'-r','Linewidth',2.0); 
hold on 
contour(Gs(2:end,:),Qs(2:end,:),dFUNDAMENTAL,[1],'-b','Linewidth',2.0); 
contour(Gs(2:end,:),Qs(2:end,:),dSECOND_HARMONIC,[1],'-g','Linewidth',2.0); 
%%%% get curves 
qOfg_threshold = 0*g0*nan;
qOfg_fundamental= 0*g0*nan;
qOfg_harmonic =0*g0*nan;

for m = 1:length(g0)
    idx_TH = find(dTHRESHOLD(:,m));
    if ~isempty(idx_TH)
    qOfg_threshold(m) =  q0(idx_TH(1)); 
    end
    idx_FUND = find(dFUNDAMENTAL(:,m));
    if ~isempty(idx_FUND)
    qOfg_fundamental(m) =  q0(idx_FUND(1)); 
    end
    idx_SEC = find(dSECOND_HARMONIC(:,m));
    if ~isempty(idx_SEC)
    qOfg_harmonic(m) =  q0(idx_SEC(1)); 
    end
end


%% do the interpolation to plot nicer stuff... 
% close
g0_2 = 1/gamma*linspace(g0(1),g0(end),1000); 

qOfg_threshold_2 = interp1(g0/gamma,qOfg_threshold,g0_2);
qOfg_fundamental_2 = interp1(g0/gamma,qOfg_fundamental,g0_2);
qOfg_harmonic_2 = interp1(g0/gamma,qOfg_harmonic,g0_2);
% the area under which fundamental mode lockign is permitted; 
FML_DIFF = qOfg_fundamental_2 - qOfg_threshold_2; 
fml_idx = find(~isnan(FML_DIFF));
FML_DIFF = FML_DIFF(fml_idx);
FML_AREA = trapz(FML_DIFF)*(g0_2(2)-g0_2(1));

HARM_DIFF = qOfg_harmonic_2 - qOfg_threshold_2;


figsize = [0,0,0.4,0.3]
figure('units','normalized','position',figsize)

plot(g0_2,qOfg_threshold_2,'-k','Linewidth',2.0);
hold on
plot(g0_2,qOfg_fundamental_2,'-b','Linewidth',2.0);
plot(g0_2,qOfg_harmonic_2,'-r','Linewidth',2.0);
legend('Threshold','Fundamental','Second Harmonic'); 
xlabel('g_0'); 
ylabel('q_0'); 

xlim([g0_2(1),g0_2(end)]); 
ylim([q0(end),q0(1)]);

f = gcf; 
for c =1:length(f.Children)
    Child = f.Children(c); 
    Child.FontName = 'Arial';
    Child.FontSize = 16; 
end










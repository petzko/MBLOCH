
INJ = settings.INJ; ULL = settings.ULL; LLL = settings.LLL; DEPOP = settings.DEPOP;

W = zeros(4);
for i = 1:4
    W(INJ,i) = settings.W_inj(i); W(ULL,i) = settings.W_ull(i);
    W(LLL,i) = settings.W_lll(i); W(DEPOP,i) = settings.W_dep(i);
end

%%%% NOTICE THE CHANGE IN RATES HERE <- we do not include the backscattering rates w13 and w12!!!
W31 = W(ULL,DEPOP) + W(ULL,INJ); W13 = W(INJ,ULL);
W21 = W(LLL,DEPOP) + W(LLL,INJ); W12 = W(INJ,LLL);
W23 = W(LLL,ULL);                W32 = W(ULL,LLL);


G1 = W13 + W12;  G3 = W31 + W32; G2 = W21 + W23;

display('no RT:');
ntot = 1;
syms n1 n2 n3
[sol1, sol3,sol2] = solve([-G1*n1 + W31*n3 + W21*n2 == 0, W13*n1-G3*n3+W23*n2 == 0, n1+n3+n2 == ntot ], [n1, n3,n2]);
double(sol1) , double(sol3) ,double(sol2)

display(['inversion: dN32 = ntot x ' num2str(double(sol3-sol2))]);
display(['resonance: dN13 = ntot x ' num2str(double(sol1-sol3))]);

%%%% rt case:

if(settings.deph>0)
    gamma_31 = 0.5*(G1+G3)+1/settings.Tdeph; %% dephsing of the resonant tunneling transition
    gamma_23 = 0.5*(G2+G3)+1/settings.Tdeph; % dephasing of the optical transision...
    gamma_21 = 0.5*(G2+G1)+1/settings.Tdeph; % dephasing of the latest transition
else
    gamma_31 = 0.5*(G1+G3); %% dephsing of the resonant tunneling transition
    gamma_23 = 0.5*(G2+G3); % dephasing of the optical transision...
    gamma_21 = 0.5*(G2+G1); % dephasing of the latest transition
end
tch = settings.tch ; lch = settings.lch;
hbar = Constants('hbar',{'time',tch})/Constants('q0');

E1 = settings.E1/hbar;
E3 = settings.E3/hbar;
E2 = settings.E2/hbar; % in rad/ps
% O13 = settings.O13/hbar; % in rad/ps

E13 = E1-E3; %rad/ps; 1->3 traisition freq
E12 = E1-E2; %rad/ps; 1->2 transition freq
E32 = E3-E2; %rad/ps  3->2 transition freq (optical transition)
% central frequency and wave number!
E0 =  (E32+E12)/2; % central OPTICAL frequency.

dE31 = 1i*E13 - gamma_31; %
dE23 = -1i*(E0 - E32) - gamma_23; %
dE21 = -1i*(E0 - E12)- gamma_21; %

n = settings.n;  c = Constants('c',{'time',tch},{'length',lch})/n;

%vacuum permitivity
Ncarriers = settings.dN*(100^3)*settings.Ld/settings.Lp; % cm^-3 --> m^-3; carrier density
Overlap = settings.Overlap;  % overlap factor -> dimensionless
k0 = (E0/c)/lch;
trace_rho = ((k0*Ncarriers*Overlap*((settings.zUL*1E-9*Constants('q0'))^2))/( Constants('eps0')*n^2*Constants('hbar')))/(1/(lch*tch));


Length = 5*1E-3;
Width = 20*1E-6;%m
Lp = 54.8*1E-9;

Anticrossings = linspace(0,5E-3,50);
N1 = []; N2 = []; N3 = []; coupling = []; J = [];
ctr  =1;

for AC = Anticrossings
    
    ctr
    O13 = AC/hbar; % in rad/ps
    syms r110 r330 r220 r310 n23p n23m n21p n21m
    [sol1,sol3,sol2,sol4,sol5,sol6,sol7,sol8] = ...
        solve( ...
        [1i*O13.*(conj(r310)-r310) + r330*W31 + r220*W21 - G1*r110 == 0 , ...
        1i*O13.*(r310 - conj(r310))+ r110*W13 + r220*W23 - G3*r330 ==0 ,...
        r110 + r330+ r220 == trace_rho, ...
        dE31*r310 - 1i*O13*(r110-r330) ==0 ,...
        dE23*n23p  + 1i*O13*n21p ==0 ,...
        dE23*n23m  + 1i*O13*n21m ==0 ,...
        dE21*n21p  + 1i*O13*n23p ==0 ,...
        dE21*n21m  + 1i*O13*n23m == 0],[r110,r330,r220,r310,n23p,n23m,n21p,n21m]);
    %         double(sol1) , double(sol3) ,double(sol2) , double(sol4)
    %         display('RT case:');
    %         display(['inversion: dN32 = ntot x ' num2str(double(sol3-sol2))]);
    %         display(['resonance: dN13 = ntot x ' num2str(double(sol1-sol3))]);
    
    N1(ctr) = double(sol1); N3(ctr) = double(sol3); N2(ctr) = double(sol2); coupling(ctr) = double(sol4);
    Jm2 = Constants('q0')*Lp*Ncarriers/trace_rho*1E12*(double(sol1)*W(INJ,DEPOP) + double(sol2)*W(LLL,DEPOP) + double(sol3)*W(ULL,DEPOP));
    J(ctr) = Jm2/(1E4);
    ctr =ctr +1;
end
% close all;
clear;
ratesfile = 'fitted_data.mat';
simDataFile = 'OPTICAMAIN.sim';
scenarioFile = 'IVL.set';


%parse all input files and load the scatterin rates file !
settings = parseSimParams(scenarioFile);
settings = parseSimData(simDataFile,settings);
settings.N =1;
load(ratesfile);

Biases = [8.0:.05:11.0]/10;
Js = 0*Biases;
global dat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZE DEVICE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%
%phase velocity inside the medium ( in mm per picosecond ... )RT time and
%frequency
dat.c = Constants('c',{'time',settings.tch},{'length',settings.lch})/settings.nTHz; dat.T_R = 2*settings.Ltot/dat.c; dat.f_R = 1/dat.T_R;
dat.dtype = 'double';
%%%%dipole mtx elements (in Cnm)
% hbar in eV-ps
dat.hbar = Constants('hbar',{'time',settings.tch})/Constants('q0');
dat.INJ = 1; dat.ULL = 2;dat.LLL = 3; dat.RES = 4; dat.DEPOP = 5;
%cavity loss l_0 in (cm^-1) --> l_0*100 in (m^-1) --> 1 mm^-1
dat.l_0 = settings.loss*100/(1/settings.lch);

dat.Ncarriers = settings.dN*(100^3)*settings.Ld/settings.Lp; % cm^-3 --> m^-3; carrier density
dat.NLVLS = 5;

dat.W = zeros(settings.N,dat.NLVLS,dat.NLVLS,dat.dtype);
dat.G = zeros(settings.N,dat.NLVLS,dat.dtype);
dat.dipR = ones(settings.N,1,dat.dtype);

RB  = 9.0/10*ones(settings.N,1);
dat.W(:,INJ,ULL) = W_fit{INJ,ULL}(RB);
dat.W(:,INJ,LLL) = W_fit{INJ,LLL}(RB);
dat.W(:,INJ,RES) = W_fit{INJ,RES}(RB);

dat.W(:,ULL,INJ) = W_fit{ULL,INJ}(RB);
dat.W(:,ULL,LLL) = W_fit{ULL,LLL}(RB);
dat.W(:,ULL,RES) = W_fit{ULL,RES}(RB);
dat.W(:,ULL,DEPOP) = W_fit{ULL,DEPOP}(RB);
%
dat.W(:,LLL,INJ) = W_fit{LLL,INJ}(RB);
dat.W(:,LLL,ULL) = W_fit{LLL,ULL}(RB);
dat.W(:,LLL,RES) = W_fit{LLL,RES}(RB);
dat.W(:,LLL,DEPOP) = W_fit{LLL,DEPOP}(RB);

dat.W(:,RES,INJ) = W_fit{RES,INJ}(RB);
dat.W(:,RES,ULL) = W_fit{RES,ULL}(RB);
dat.W(:,RES,LLL) = W_fit{RES,LLL}(RB);
dat.W(:,RES,DEPOP) = W_fit{RES,DEPOP}(RB);

dat.G(:,INJ) = dat.W(:,INJ,ULL)+ dat.W(:,INJ,LLL)+dat.W(:,INJ,RES);
dat.G(:,ULL) = dat.W(:,ULL,INJ)+dat.W(:,ULL,LLL)+dat.W(:,ULL,DEPOP) + dat.W(:,ULL,RES);
dat.G(:,LLL) =  dat.W(:,LLL,INJ)+dat.W(:,LLL,ULL)+dat.W(:,LLL,DEPOP) + dat.W(:,LLL,RES) ;
dat.G(:,RES) =  dat.W(:,RES,INJ)+dat.W(:,RES,ULL)+dat.W(:,RES,LLL) +dat.W(:,RES,DEPOP) ;
% % % % Added Pure Dephasing
if(settings.deph>0)
    dat.gamma_13 = 0.5*(dat.G(:,dat.INJ)+dat.G(:,dat.ULL))+1/settings.Tdeph_1; %% dephsing of the resonant tunneling transition
    dat.gamma_32 = 0.5*(dat.G(:,dat.ULL)+dat.G(:,dat.LLL))+1/settings.Tdeph_2; % dephasing of the optical transision...
    dat.gamma_12 = 0.5*(dat.G(:,dat.INJ)+dat.G(:,dat.LLL))+1/settings.Tdeph_3; % dephasing of the latest transition
else
    dat.gamma_13 = 0.5*(dat.G(:,dat.INJ)+dat.G(:,dat.ULL)); %% dephsing of the resonant tunneling transition
    dat.gamma_32 = 0.5*(dat.G(:,dat.ULL)+dat.G(:,dat.LLL)); % dephasing of the optical transision...
    dat.gamma_12 = 0.5*(dat.G(:,dat.INJ)+dat.G(:,dat.LLL)); % dephasing of the latest transition
end

  X = zeros(15,1);
%     X(1) = 1/3;
%     X(2) = 1/3;
%     X(3) = 1/3;
%     X(4) = 0;
%     X(5) = 0;
%     X(6) = 0;
    
    f =0.001; 
%     
%     X(7) = f*(rand(1)+1i*rand(1));
%     X(8) = f*(rand(1)+1i*rand(1));
%     X(9) = f*(rand(1)+1i*rand(1));
%     X(10) = f*(rand(1)+1i*rand(1));
%     X(11) = f*(rand(1)+1i*rand(1));
%     X(12) = f*(rand(1)+1i*rand(1));
%     X(13) = f*(rand(1)+1i*rand(1));
%     X(14) = f*(rand(1)+1i*rand(1));
%     X(15) = f*(rand(1)+1i*rand(1));

    
for ctr_ = 1:length(Biases)
    bias = Biases(ctr_);
    
    dat.v_TL = bias;
    dat.diffusion = 0 ;
    %central field frequency.
    dat.E0 =  (E_fit{dat.ULL}(bias)-E_fit{dat.LLL}(bias))/dat.hbar;
    dat.zUL = zUL_fit(bias);
    dat.trace_rho = ((dat.E0*1E12*dat.Ncarriers*settings.Overlap*((dat.zUL*1E-9*Constants('q0'))^2))/(Constants('eps0')*settings.nTHz*Constants('c')*Constants('hbar')))/(1/(settings.lch*settings.tch));
    
    dat = interpParams(settings,dat,W_fit,E_fit,AC_fit,zUL_fit);
    
    
    dat.losses = -dat.l_0;
    dat.factor = -1i*dat.trace_rho;
    
    
    X(1) = 1/4;
    X(2) = 1/4;
    X(3) = 1/4;
    X(4) = 1/4;
    X(5) = 0;
    X(6) = 0;
    
    f = 1E-5; 
    X(7) = f*rand(1)*((1)+1i*(1));
    X(8) = f*rand(1)*((1)-1i*(1));
    X(9) = f*rand(1)*((1)-1i*(1));
    X(10) = f*rand(1)*((1)-1i*(1));
    X(11) = f*rand(1)*((1)-1i*(1));
    X(12) = f*rand(1)*((1)-1i*(1));
    X(13) = f*rand(1)*(1-1i*(1));
    X(14) = 0.001+1i*0.0001;
    X(15) = 0.001+1i*0.0001;
    
    
    % now solve the DM equations
    
    W = squeeze(dat.W);
    options = optimoptions('fsolve','Display','none');
    [sol,val] = fsolve(@FP_DM_RHS,X,options);
    
    % X(1) = r110,
    % X(2) = r330;
    % X(3) = r220;
    % X(4) = r11p;
    % X(5) = r33p;
    % X(6) = r22p;
    % X(7) = rRES;
    % X(8) = r130;
    % X(9) = r13p;
    % X(10) = r13m;
    % X(11) = n32p;
    % X(12) = n32m;
    % X(13) = n12p;
    % X(14) = n12m;
    
    r13 =  sol(7) +(sol(8) +sol(9));
    rates1 = 1i*dat.O13.*(r13-conj(r13));
    
    r11 = sol(1) + sol(4)+conj(sol(4));
    r33 = sol(2) + sol(5)+conj(sol(5));
    r22 = sol(3) + sol(6)+conj(sol(6));
    rRES = 1-sol(1)-sol(2)-sol(3);
    rates2 =  r11.*W(INJ,DEPOP) + r33.*W(ULL,DEPOP)+r22.*W(LLL,DEPOP) +rRES.*W(RES,DEPOP);
    
    sum = sol(1)+sol(2)+sol(3)
    
    Js1(ctr_) =-(Constants('q0')*(settings.Lp*1E-9)*dat.Ncarriers*1E12*rates1)/1E6; %in A/cm^2
    Js2(ctr_) =(Constants('q0')*(settings.Lp*1E-9)*dat.Ncarriers*1E12*rates2)/1E6; %in A/cm^2
    L(ctr_) = abs(sol(end)).^2+abs(sol(end-1)); 
    % end of solving the DM equaitons
    
%     X = sol;
    
    
end

% figure;
hold on;
plotyy(Biases(3:end),Js1(3:end)*3,Biases(3:end),L(3:end));
% plot(Biases,real(Js2));

function [ savename ] = runsim_MB_TL_3lvl( scenarioFile,simDataFile,ratesfile,varargin )


close all;

lvar = length(varargin);

%handle user input!
if(length(varargin)>0)
    if mod(lvar,2) ~=0
        display('incorrect number of input arguments. Aborting!');
        return;
    end
    for i = 1:2:lvar
        name = varargin{i}; val  = varargin{i+1};
        name = strtrim(lower(name)); %change to lower case!
        if strcmp(name,'init')
            init = val;
        else if strcmp(name,'workspace')
                %regexp to load all variables except init!
                load(val,'-regexp','^(?!.*init.*).*$');
                
                
            else
                display( ['Unknown input option name: ' name '. Aborting!' ]);
                return;
            end
        end
    end
else
    init = 1;
end

%parse all input files and load the scatterin rates file !
settings = parseSimParams(scenarioFile);
settings = parseSimData(simDataFile,settings);
load(ratesfile);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEVICE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%
%phase velocity inside the medium ( in mm per picosecond ... )RT time and
%frequency 
dat.c = Constants('c',{'time',settings.tch},{'length',settings.lch})/settings.nTHz; dat.T_R = 2*settings.Ltot/dat.c; dat.f_R = 1/dat.T_R;

%%%%dipole mtx elements (in Cnm)
% hbar in eV-ps
dat.hbar = Constants('hbar',{'time',settings.tch})/Constants('q0');

dat.INJ = 2; 
dat.DEPOP = 6;
dat.ULL = 3;
dat.LLL = 4;


dat.Ncarriers = settings.dN*(100^3)*settings.Ld/settings.Lp; % cm^-3 --> m^-3; carrier density

%calculate the normalization constant, i.e. the number on which we divide
%all quantities ( except the electric field envelope) to simplify
%our initial system. it is important as it determines the initial value of
%the overall electron population inside the system-> a quantity that shall
%be perserved throughout the whole simulaiton ! !

dat.E0 =  3.8666*2*pi;%(E32+E12)/2; % central OPTICAL frequency.

dat.trace_rho = ((dat.E0*1E12*dat.Ncarriers*settings.Overlap*((settings.zUL*1E-9*Constants('q0'))^2))/(Constants('eps0')*settings.nTHz*Constants('c')*Constants('hbar')))/(1/(settings.lch*settings.tch));
%cavity loss l_0 in (cm^-1) --> l_0*100 in (m^-1) --> 1 mm^-1
dat.l_0 = settings.loss*100/(1/settings.lch);

%%%%%%%%%%%%%%%%% simuation parameters %%%%%%%%%%%%%%%%%

%grid size in x direction
dat.x = linspace(0,settings.Ltot,settings.N)';
dat.dx = dat.x(2) - dat.x(1); dat.dt = dat.dx/dat.c;

dat.diffusion = 4*dat.E0/dat.c^2*settings.D*10^2/(1/settings.tch);

%%%% specify some of the mainloop control parameters %%%%
iter_per_rt = round(dat.T_R/dat.dt);
checkptIter = iter_per_rt*100; %make a checkpoint every 100RTs.

tEnd = settings.simRT*dat.T_R; % end time in tps

plotCtr = settings.plotCtr; %% set on how many iterations should the program plot the intensity

%simulation info storage arrays -> preallocate
recordingduration = settings.recordRT*dat.T_R; % how many ps should we record the pulse for

if(init > 0)
    
    
    clc;
    
    W = settings.Wmtx;
    dat.NLVLS = sqrt(length(W)); % nr of levels to consider !
    %reshape into a matrix
    dat.W = reshape(W,dat.NLVLS,dat.NLVLS).';
    
    
    
    dat.G = zeros(dat.NLVLS,1);
    dat.i0 =  17%A/mm
    
    dat = makeMaxwellVars(settings,dat);
    dat = makeBlochVars(settings,dat);  
    dat = makeTransLineVars(settings,dat); 
    t =1;
    dat.t = dat.dt; 
    
    idx = settings.N; % the index of the predefined point we sample the feild at!
    iter_ctr = 0;  ctr = 1;
    
    record_U= 1;
    record_V = 1; 

    record_r110 = 1; 
    record_r330 = 1; 
    record_r220 = 1;
    

    record_v_TL =1; 
    record_i_TL =1;
    record_J_TL =1;
    
    
      
    %convert the curr density in A/m^2 into A/mm^2 by dividing by 1E6;
    
    
end

iterperrecord = 1;
recordingiter  = round(recordingduration/iterperrecord/dat.dt);
padsize = recordingiter-length(record_U);

%preallocate memory for optical field, curr density voltage wave, etc.
%strorage...
record_U= padarray(record_U,padsize,'post'); record_V = padarray(record_V,padsize,'post');

record_v_TL = padarray(record_v_TL,padsize,'post'); 
record_i_TL = padarray(record_i_TL,padsize,'post');
record_J_TL = padarray(record_J_TL,padsize,'post');

%store population info
record_r110 = padarray(record_r110,padsize,'post');    
record_r330 = padarray(record_r330,padsize,'post');
record_r220 = padarray(record_r220,padsize,'post');    

info.settings = settings;
info.cavity = 'FP';
info.Ltot = settings.Ltot;
info.N = settings.N;
info.SIMTYPE = 'WITHOUT DISP. COMP';

dat.G(dat.INJ) = dat.W(dat.INJ,dat.ULL)+ dat.W(dat.INJ,dat.LLL);
dat.G(dat.ULL) = dat.W(dat.ULL,dat.INJ)+dat.W(dat.ULL,dat.LLL)+dat.W(dat.ULL,dat.DEPOP);
dat.G(dat.LLL) =  dat.W(dat.LLL,dat.INJ)+dat.W(dat.LLL,dat.ULL)+dat.W(dat.LLL,dat.DEPOP);

   % % % % Added Pure Dephasing
if(settings.deph>0)
    dat.gamma_13 = 0.5*(dat.G(dat.INJ)+dat.G(dat.ULL))+1/settings.Tdeph_1; %% dephsing of the resonant tunneling transition
    dat.gamma_32 = 0.5*(dat.G(dat.LLL)+dat.G(dat.ULL))+1/settings.Tdeph_2; % dephasing of the optical transision...
    dat.gamma_12 = 0.5*(dat.G(dat.LLL)+dat.G(dat.INJ))+1/settings.Tdeph_3; % dephasing of the latest transition
else
    dat.gamma_13 = 0.5*(dat.G(dat.INJ)+dat.G(dat.ULL)); %% dephsing of the resonant tunneling transition
    dat.gamma_32 = 0.5*(dat.G(dat.LLL)+dat.G(dat.ULL)); % dephasing of the optical transision...
    dat.gamma_12 = 0.5*(dat.G(dat.LLL)+dat.G(dat.INJ)); % dephasing of the latest transition
end  

while(dat.t< tEnd)
    
    %%%%% correctly setup the energies, the anticrossings and the scattering rates...
    %obtain the energies of the core levels for the simulation
    E1 = E_fit{dat.INJ}(dat.v_TL)/dat.hbar;
    E3 = E_fit{dat.ULL}(dat.v_TL)/dat.hbar;
    E2 = E_fit{dat.LLL}(dat.v_TL)/dat.hbar; % in rad/ps
    
%     E1 = E3;
    
    %%%% energies: E1 > E3 > E2 Optical transition is 3->2 TUNNELING transition is 1->3
    dat.O13 = AC_fit{2}(dat.v_TL)/dat.hbar; % in rad/ps
    dat.E13 = E1-E3; %rad/ps; 1->3 traisition freq
    dat.E12 = E1-E2; %rad/ps; 1->2 transition freq
    dat.E32 = E3-E2; %rad/ps  3->2 transition freq (optical transition)
    
    dat.dE13 = -1i*dat.E13 - dat.gamma_13; %
    dat.dE32 = +1i*(dat.E0 - dat.E32) - dat.gamma_32; %
    dat.dE12 = +1i*(dat.E0 - dat.E12)- dat.gamma_12; %
        
    %%%%% end of setting up the TL PARAMS %%%%%
    
    
    if(mod(iter_ctr+1,checkptIter) == 0 )
        checkptname = ['CHCKPT_' settings.name '_' settings.scenario '_N_TRANSMISSION_LINE_' num2str(settings.N) '_FP'];
        save(checkptname);
    end
    
    %%plot some of the results if neeed ariseth :D
    if(mod(iter_ctr,100) == 0)
        clc;
        info.iter_ctr = iter_ctr;
        info.RT = dat.t/dat.T_R;
        intensity = dat.U.*conj(dat.U) + dat.V.*conj(dat.V) ;
        info.maxInt  =  max(intensity);
        printINFO(info);
        
        subplot(3,1,1)
        plot(dat.x,[real(dat.U),real(dat.V)]);
        subplot(3,1,2)
        plotyy(dat.x,[dat.v_TL*10],dat.x,dat.i_TL);
        subplot(3,1,3)
        plot(dat.x,dat.J_TL);

        getframe;
    end
    %%%% obtain the field, field intensity and the total population at position "idx" ...
    
    if((dat.t >= tEnd - recordingduration) && mod(iter_ctr,iterperrecord) == 0)
        
        %store fields info
        record_U(ctr)= dat.U(idx);  record_V(ctr)= dat.V(idx);
        record_v_TL(ctr) = dat.v_TL(idx);
        record_i_TL(ctr) = dat.i_TL(idx);
        record_J_TL(ctr) = dat.J_TL(idx);
    
        record_r110(ctr) = dat.r110(idx);
        record_r220(ctr) = dat.r220(idx);
        record_r330(ctr) = dat.r330(idx);
        
        ctr = ctr+1;
    end
    
    dat = stepBloch(settings,dat);
    dat = stepWave(settings,dat);
    dat = stepTransLine(settings,dat);
    dat = updateBloch(settings,dat);
    
    
    dat.t = dat.t+dat.dt;
    iter_ctr = iter_ctr + 1;
    
end
savename = [settings.name '_' settings.scenario '_N_' num2str(settings.N) '_FP_' num2str(settings.simRT) ];
save(savename);
end

simDataFile = 'OPTICAMAIN.sim';
scenarioFile = '12-1-off.set';

%parse all input files and load the scatterin rates file !
settings = parseSimParams(scenarioFile);
settings = parseSimData(simDataFile,settings);
init = 1;
if(init > 0)
    
    clc
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZE DEVICE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%
    %phase velocity inside the medium ( in mm per picosecond ... )RT time and
    %frequency
    dat.c = Constants('c',{'time',settings.tch},{'length',settings.lch})/settings.nTHz; dat.T_R = 2*settings.Ltot/dat.c; dat.f_R = 1/dat.T_R;
    dat.dtype = 'double';
    %%%%dipole mtx elements (in Cnm)
    % hbar in eV-ps
    dat.hbar = Constants('hbar',{'time',settings.tch})/Constants('q0');
    dat.INJ = 1; dat.ULL = 2;dat.LLL = 3; dat.RES = 4; dat.DEPOP = 5;
    
    dat.NLVLS = 5;
    dat.G = zeros(dat.NLVLS,1,dat.dtype);
    dat.dipR = ones(1,1,dat.dtype);
    
    
    dat.HTB = zeros(dat.NLVLS,dat.NLVLS,dat.dtype);
    dat.W = zeros(dat.NLVLS,dat.NLVLS,dat.dtype);
    for i = 1:dat.NLVLS
        for j = 1:dat.NLVLS
            dat.HTB(i,j) = settings.HTB(j+(i-1)*dat.NLVLS);
            dat.W(i,j) = settings.Wmtx(j+(i-1)*dat.NLVLS);
        end
    end
    
    %central field frequency.
    dat.E0 =  (dat.HTB(dat.ULL,dat.ULL)-dat.HTB(dat.LLL,dat.LLL))/dat.hbar;
    
    
    %%%%%%%%%%%%%%%%%
    dat.G(dat.INJ) = dat.W(dat.INJ,dat.ULL)+ dat.W(dat.INJ,dat.LLL)+dat.W(dat.INJ,dat.RES);
    dat.G(dat.ULL) = dat.W(dat.ULL,dat.INJ)+dat.W(dat.ULL,dat.LLL)+dat.W(dat.ULL,dat.DEPOP) + dat.W(dat.ULL,dat.RES);
    dat.G(dat.LLL) =  dat.W(dat.LLL,dat.INJ)+dat.W(dat.LLL,dat.ULL)+dat.W(dat.LLL,dat.DEPOP) + dat.W(dat.LLL,dat.RES) ;
    dat.G(dat.RES) =  dat.W(dat.RES,dat.INJ)+dat.W(dat.RES,dat.ULL)+dat.W(dat.RES,dat.LLL) +dat.W(dat.RES,dat.DEPOP) ;
    
    % % % Added Pure Dephasing
    if(settings.deph>0)
        dat.gamma_13 = 0.5*(dat.G(dat.INJ)+dat.G(dat.ULL))+1/settings.Tdeph_1; %% dephsing of the resonant tunneling transition
        dat.gamma_32 = 0.5*(dat.G(dat.ULL)+dat.G(dat.LLL))+1/settings.Tdeph_2; % dephasing of the optical transision...
        dat.gamma_12 = 0.5*(dat.G(dat.INJ)+dat.G(dat.LLL))+1/settings.Tdeph_3; % dephasing of the latest transition
    else
        dat.gamma_13 = 0.5*(dat.G(dat.INJ)+dat.G(dat.ULL)); %% dephsing of the resonant tunneling transition
        dat.gamma_32 = 0.5*(dat.G(dat.ULL)+dat.G(dat.LLL)); % dephasing of the optical transision...
        dat.gamma_12 = 0.5*(dat.G(dat.INJ)+dat.G(dat.LLL)); % dephasing of the latest transition
    end
    
    %%%%% correctly setup the energies, the anticrossings and the scattering rates...
    %obtain the energies of the core levels for the simulation
    E1 = dat.HTB(dat.INJ,dat.INJ)/dat.hbar;
    E3 = dat.HTB(dat.ULL,dat.ULL)/dat.hbar;
    E2 = dat.HTB(dat.LLL,dat.LLL)/dat.hbar; % in rad/ps
    
    
    
    dat.E13 = E1-E3; %rad/ps; 1->3 traisition freq
    dat.E12 = E1-E2; %rad/ps; 1->2 transition freq
    dat.E32 = E3-E2; %rad/ps  3->2 transition freq (optical transition)
    
    dat.dE13 = -1i*dat.E13 - dat.gamma_13; %
    dat.dE32 = +1i*(dat.E0 - dat.E32) - dat.gamma_32; %
    dat.dE12 = +1i*(dat.E0 - dat.E12)- dat.gamma_12; %
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %cavity loss l_0 in (cm^-1) --> l_0*100 in (m^-1) --> 1 mm^-1
    dat.l_0 = settings.loss*100/(1/settings.lch);
    
    %%%%%%%%%%%%%%%%% simuation parameters %%%%%%%%%%%%%%%%%
    %grid size in x direction
    dat.dx = settings.Ltot/(settings.N-1);
    dat.x = [0:settings.N-1].'*dat.dx;
    
    if strcmp(dat.dtype,'single')
        dat.x = single(dat.x);
    end
    
    dat.dt = dat.dx/dat.c;
    dat.diffusion = 4*dat.E0/dat.c^2*settings.D*10^2/(1/settings.tch);
    
    dat.zUL = settings.zUL;
    dat.Ncarriers = 5./3.*settings.dN*(100^3)*settings.Ld/settings.Lp; % cm^-3 --> m^-3; carrier density
    dat.trace_rho = ((dat.E0*1E12*dat.Ncarriers*settings.Overlap*((dat.zUL*1E-9*Constants('q0'))^2))/(Constants('eps0')*settings.nTHz*Constants('c')*Constants('hbar')))/(1/(settings.lch*settings.tch));
    
    
    
    
    dat.t = dat.dt;
    idx = 100; % the index of the predefined point we sample the feild at!
    dat.iter_ctr = 0;  ctr = 1;
    
    record_U= 1;
    record_V = 1;
    
    record_r110 = 1;
    record_r330 = 1;
    record_r220 = 1;
    record_rRES = 1;
    
    
    record_popsum = 1;
    record_v_TL =1;
    record_i_TL = 1;
    record_J_TL =1;
    
    dat = makeMaxwellVars(settings,dat);
    dat = makeBlochVars(settings,dat);
    
end

%%%% specify some of the mainloop control parameters %%%%
iter_per_rt = round(dat.T_R/dat.dt);
checkptIter = iter_per_rt*100; %make a checkpoint every 100RTs.

N_t  = 10000; % iter_per_rt*settings.simRT;

plotCtr = settings.plotCtr; %% set on how many iterations should the program plot the intensity
interpCtr = 100; %set how often to interpolate the energies, scattering rates, dipole elements etc.

%simulation info storage arrays -> preallocate
recordingduration = N_t; % how many ps should we record the pulse for

iterperrecord = 1;
recordingiter  = round(recordingduration/iterperrecord/dat.dt);
padsize = double(recordingiter-length(record_U));

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
record_rRES = padarray(record_rRES,padsize,'post');

record_popsum = padarray(record_popsum,padsize,'post');

if(strcmp(dat.dtype,'single'))
    
    record_U= single(record_U); record_V = single(record_V);
    
    %store population info
    record_r110 = single(record_r110);
    record_r330 = single(record_r330);
    record_r220 = single(record_r220);
    record_rRES = single(record_rRES);
    
    record_popsum = single(record_popsum);
    
end

info.settings = settings;
info.cavity = 'FP-OPTICA-NOTL';
info.Ltot = settings.Ltot;
info.N = settings.N;
info.SIMTYPE = ['OPTICA IN FP CAVITY @ 9.8 kV/cm'];
tic;
totaldur =0; 
while(dat.iter_ctr< N_t)
    
  
    %%%%% end of setting up the TL PARAMS %%%%%
    
    if(mod(dat.iter_ctr+1,checkptIter) == 0 )
        checkptname = ['CHCKPT_' settings.name '_' settings.scenario '_N_TRANSMISSION_LINE_' num2str(settings.N) '_FP'];
        save(checkptname);
    end
    
    %%plot some of the results if neeed ariseth :D
    if(mod(dat.iter_ctr,1000) == 0)
        te = toc;
        totaldur = totaldur+te;
        tic;
        clc;
        info.iter_ctr = dat.iter_ctr;
        info.RT = dat.t/dat.T_R;
        intensity = dat.U.*conj(dat.U) + dat.V.*conj(dat.V) ;
        info.maxInt  =  max(intensity);
        printINFO(info);
        
%         subplot(2,1,1);
%         plot(dat.x,[real(dat.U),real(dat.V)]);
%         
%         title(info.SIMTYPE);
%         subplot(2,1,2)
%         %plots the populations and the current density
%         plot(dat.x,[dat.r110,dat.r330,dat.r220,dat.rRES]);
%         getframe;
        display(['calculation time (s): ' num2str(te)]);

        
    end
    
    %%%% obtain the field, field intensity and the total population at position "idx" ...
    
        
    %store fields info
    record_U(ctr)= dat.U(idx);  
    record_V(ctr)= dat.V(idx);
    record_r110(ctr) = dat.r110(idx);
    record_r220(ctr) = dat.r220(idx);
    record_r330(ctr) = dat.r330(idx);
    record_rRES(ctr) = dat.rRES(idx);

    popsum = dat.r110(idx)+dat.r220(idx)+dat.r330(idx) + dat.rRES(idx);
    record_popsum(ctr) =  popsum;
    ctr = ctr+1;
    
    
    dat = stepBloch(settings,dat);
    dat = stepWave(settings,dat);
    dat = updateBloch(settings,dat);
    
    dat.t = dat.t+dat.dt;
    dat.iter_ctr = dat.iter_ctr + 1;
    
end

display(['Simulation time: ' num2str(dat.t) ' (ps)'])
display(['-> calculation real time: (sec) ' num2str(totaldur) ])

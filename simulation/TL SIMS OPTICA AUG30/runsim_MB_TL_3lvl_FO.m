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
                load(val,'-regexp','-regexp','^(?!(init|scenarioFile|simDataFile|varargin)$).');
                
                
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
    
    %central field frequency.
    dat.E0 =  (E_fit{dat.ULL}(settings.bias/10)-E_fit{dat.LLL}(settings.bias/10))/dat.hbar;
    %cavity loss l_0 in (cm^-1) --> l_0*100 in (m^-1) --> 1 mm^-1
    dat.l_0 = settings.loss*100/(1/settings.lch);
    
    %%%%%%%%%%%%%%%%% simuation parameters %%%%%%%%%%%%%%%%%
    %grid size in x direction
    dat.x = linspace(0,settings.Ltot,settings.N)';
    if strcmp(dat.dtype,'single')
        dat.x = single(dat.x);
    end
    
    dat.dx = dat.x(2) - dat.x(1); dat.dt = dat.dx/dat.c;
    dat.diffusion = 4*dat.E0/dat.c^2*settings.D*10^2/(1/settings.tch);
    
    dat.zUL = zUL_fit(settings.bias/10);
    dat.Ncarriers = settings.dN*(100^3)*settings.Ld/settings.Lp; % cm^-3 --> m^-3; carrier density
    dat.trace_rho = ((dat.E0*1E12*dat.Ncarriers*settings.Overlap*((dat.zUL*1E-9*Constants('q0'))^2))/(Constants('eps0')*settings.nTHz*Constants('c')*Constants('hbar')))/(1/(settings.lch*settings.tch));
    
    
    dat.NLVLS = 5;
    dat.W = zeros(settings.N,dat.NLVLS,dat.NLVLS,dat.dtype);
    dat.G = zeros(settings.N,dat.NLVLS,dat.dtype);
    dat.dipR = ones(settings.N,1,dat.dtype);
    
    
    dat.t = dat.dt;
    idx = settings.N; % the index of the predefined point we sample the feild at!
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
    
    
    
    dat.v_TL = settings.bias/10*ones(settings.N,1,dat.dtype);
    dat = interpParams(settings,dat,W_fit,E_fit,AC_fit,zUL_fit);
    dat = makeMaxwellVars(settings,dat);
    dat = makeTransLineVarsx2(settings,dat);
    dat = makeBlochVars(settings,dat);
    
end

%%%% specify some of the mainloop control parameters %%%%
iter_per_rt = round(dat.T_R/dat.dt);
checkptIter = iter_per_rt*100; %make a checkpoint every 100RTs.

tEnd = settings.simRT*dat.T_R; % end time in tps
plotCtr = settings.plotCtr; %% set on how many iterations should the program plot the intensity
interpCtr = 500; %set how often to interpolate the energies, scattering rates, dipole elements etc.

%simulation info storage arrays -> preallocate
recordingduration = settings.recordRT*dat.T_R; % how many ps should we record the pulse for

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
    
    record_v_TL = single(record_v_TL);
    record_i_TL = single(record_i_TL);
    record_J_TL = single(record_J_TL);
    
    %store population info
    record_r110 = single(record_r110);
    record_r330 = single(record_r330);
    record_r220 = single(record_r220);
    record_rRES = single(record_rRES);
    
    record_popsum = single(record_popsum);
    
end

info.settings = settings;
info.cavity = 'FP-OPTICA';
info.Ltot = settings.Ltot;
info.N = settings.N;
info.SIMTYPE = ['TL with (i0,v0) = (' num2str(settings.current*10) ',' num2str(settings.bias) ') in units (A/cm,kV/cm)'];


while(dat.t< tEnd)
    
    if(mod(dat.iter_ctr+1,interpCtr) == 0)
        dat = interpParams(settings,dat,W_fit,E_fit,AC_fit,zUL_fit);
    end
    
    %%%%% end of setting up the TL PARAMS %%%%%
    
    if(mod(dat.iter_ctr+1,checkptIter) == 0 )
        checkptname = ['CHCKPT_' settings.name '_' settings.scenario '_N_TRANSMISSION_LINE_' num2str(settings.N) '_FP'];
        save(checkptname);
    end
    
    %%plot some of the results if neeed ariseth :D
    if(mod(dat.iter_ctr,100) == 0)
        clc;
        info.iter_ctr = dat.iter_ctr;
        info.RT = dat.t/dat.T_R;
        intensity = dat.U.*conj(dat.U) + dat.V.*conj(dat.V) ;
        info.maxInt  =  max(intensity);
        printINFO(info);
        display(['Vs: ' num2str(dat.Vs(dat.t))]);
        display(['v_0: ' num2str(dat.v_TL(1))]);
        display(['i-in: ' num2str(dat.i_TL(1))]);
        display(['i-out_1: ' num2str(trapz(dat.x,dat.J_TL1))]);
        display(['i-out_2: ' num2str(trapz(dat.x,dat.J_TL2))]);
        
        
        subplot(3,1,1)
        plot(dat.x,[real(dat.U),real(dat.V)]);
        subplot(3,1,2)
        ax = plotyy(dat.x,dat.v_TL*10,dat.x(1:end-1),dat.i_TL);

        title(info.SIMTYPE);
        subplot(3,1,3)
        %plots the populations and the current density
        plotyy(dat.x,[dat.r110,dat.r330,dat.r220,dat.rRES],dat.x,dat.J_TL1);
        getframe;
%         
    end
    %%%% obtain the field, field intensity and the total population at position "idx" ...
    
    if(mod(dat.iter_ctr,iterperrecord) == 0)
        
        %store fields info
        record_U(ctr)= dat.U(idx);  record_V(ctr)= dat.V(idx);
        record_v_TL(ctr) = dat.v_TL(1);
        record_i_TL(ctr) = dat.i_TL(1);
        
        record_J_TL(ctr) = dat.J_TL1(idx);
        
        record_r110(ctr) = dat.r110(idx);
        record_r220(ctr) = dat.r220(idx);
        record_r330(ctr) = dat.r330(idx);
        record_rRES(ctr) = dat.rRES(idx);
        
        popsum = dat.r110(idx)+dat.r220(idx)+dat.r330(idx) + dat.rRES(idx);
        record_popsum(ctr) =  popsum;
        ctr = ctr+1;
    end
    
    dat = stepBloch(settings,dat);
    dat = stepWave(settings,dat);
    dat = updateBloch(settings,dat);    
    dat = stepTransLinex2(settings,dat);
    
    
    dat.t = dat.t+dat.dt;
    dat.iter_ctr = dat.iter_ctr + 1;
    
end
savename = [settings.name '_' settings.scenario '_N_' num2str(settings.N) '_FP_' num2str(settings.simRT) ];
save(savename);
end

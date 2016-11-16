
close all;
simDataFile = 'SITTEST.sim'
%parse all input files and load the scatterin rates file !
sim_settings = parseSimDataSIT(simDataFile);

init = 1

dat.tp = .1;
npulse = 1; 
dat.A0 = 4/dat.tp;
dat.t0 = dat.tp*10;

if(init > 0)
    
    clc
    %length of the absorption and the gain region
    dat.Ltot =  sim_settings.Lg + sim_settings.La;
    
    %grid size of the absorption section
    dat.Na = round((sim_settings.La/(sim_settings.Lg+sim_settings.La))*sim_settings.N);
   
    
    %grid size of the gain section
    dat.Ng = sim_settings.N-dat.Na;
    
    
    % ordering: absorption section -> gain section
    dat.abs_start = 1; dat.abs_end = dat.Na;
    dat.gain_start = dat.Na +1 ; dat.gain_end = sim_settings.N;
    % ordering: gain section -> absorption section
    % dat.gain_start = 1; dat.gain_end = dat.Ng;
    % dat.abs_start = dat.Ng +1 ; dat.abs_end = settings.N;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZE DEVICE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%
    % central freq. phase velocity in characteristic units
    dat.c = Constants('c',{'time',sim_settings.tch},{'length',sim_settings.lch})/sim_settings.nTHz;
    % round trip time and frequency
    dat.T_R = dat.Ltot/dat.c; dat.f_R = 1/dat.T_R;
    % simulation datatype (double or single)
    dat.dtype = 'double';
    % hbar in eV-ps
    dat.hbar = Constants('hbar',{'time',sim_settings.tch})/Constants('q0');
    
    %grid size in x direction
    dat.x = linspace(0,dat.Ltot,sim_settings.N)';
    if strcmp(dat.dtype,'single')
        dat.x = single(dat.x);
    end
    dat.dx = dat.x(2) - dat.x(1); dat.dt = dat.dx/dat.c;
     
    dat.core = dat.x >= dat.Ltot/20 & dat.x <= 19*dat.Ltot/20;
    
    
    
    %indices of the upper and lower laser level
    dat.ULL = sim_settings.ULL ; dat.LLL = sim_settings.LLL;
    %indices of the gain and absorptopn secions
    dat.GAIN = 1; dat.ABS = 2;
    %eigenenergies of the gain hamiltonian -> transform into ang. freq.
    dat.NLVLS = 2;
    Hg = reshape(sim_settings.Hg',dat.NLVLS,dat.NLVLS);
    dat.E2g = Hg(dat.ULL,dat.ULL)/dat.hbar;
    dat.E1g = Hg(dat.LLL,dat.LLL)/dat.hbar;
    
    %eigenenergies of the gain hamiltonian -> transform into ang. freq.
    Ha = reshape(sim_settings.Ha',dat.NLVLS,dat.NLVLS);
    dat.E2a = Ha(dat.ULL,dat.ULL)/dat.hbar;
    dat.E1a = Ha(dat.LLL,dat.LLL)/dat.hbar;
    
    dat.E21g = dat.E2g-dat.E1g; %rad/ps  2->1 gain transition freq (optical transition)
    dat.E21a = dat.E2a-dat.E1a; %rad/ps; 2->1 absorber transition freq
    %central field frequency.
    dat.E0 = (dat.E21g+dat.E21a)/2;
    
    %cavity loss l_0 in (cm^-1) --> l_0*100 in (m^-1) --> 1 mm^-1
    dat.lg = sim_settings.gain_loss*100/(1/sim_settings.lch);
    dat.la = sim_settings.abs_loss*100/(1/sim_settings.lch);
    
    %gain section carrier density (doping density) and trace normalization
    %constant
    dat.Ncarriers_g = sim_settings.dN_g*(100^3)*sim_settings.Ld_g/sim_settings.Lp_g; % cm^-3 --> m^-3; carrier density
    dat.trace_rho_g = ((dat.E0*1E12*dat.Ncarriers_g*sim_settings.Overlap_g* ...
        ((sim_settings.zULg*1E-9*Constants('q0'))^2))/(Constants('eps0')*sim_settings.nTHz*Constants('c')*Constants('hbar')))/(1/(sim_settings.lch*sim_settings.tch));
    
    %absorption section carrier density and trace normalization constant
    dat.Ncarriers_a = sim_settings.dN_a*(100^3)*sim_settings.Ld_a/sim_settings.Lp_a; % cm^-3 --> m^-3; carrier density
    dat.trace_rho_a = ((dat.E0*1E12*dat.Ncarriers_a*sim_settings.Overlap_a* ...
        ((sim_settings.zULa*1E-9*Constants('q0'))^2))/(Constants('eps0')*sim_settings.nTHz*Constants('c')*Constants('hbar')))/(1/(sim_settings.lch*sim_settings.tch));
    
    
    dat.NLVLS = 2;

    % dephasing of the latest transition
    %implementation # 2 -> rostislav model ! 
    dat.gamma_21_g = 1./sim_settings.T2_g;
    dat.gamma_21_a = 1./sim_settings.T2_a;
    
    dat.dE21 = ones(sim_settings.N,1,dat.dtype);
    dat.dE21(dat.gain_start:dat.gain_end) = +1i*(dat.E0 - dat.E21g) - ...
        dat.gamma_21_g; %
    dat.dE21(dat.abs_start:dat.abs_end) = +1i*(dat.E0 - dat.E21a) - ...
        dat.gamma_21_a; %
    
    %the varying dipole ratio.
    dat.dipR = ones(sim_settings.N,1,dat.dtype);
    dat.dipR(dat.abs_start:dat.abs_end) = sim_settings.zULg/sim_settings.zULa;
    dat.dipR_inv = 1./dat.dipR;
    
    dat = makeMaxwellVars(sim_settings,dat);
    dat = makeBlochVars(sim_settings,dat);
    
    
    dat.t = dat.dt;
    % the index of the predefined point iun the
    % gain and absorption media where to sample the feild at!
    idx_g = dat.gain_start;
    idx_a = dat.abs_start;
    idx = 1; 
    
    
    iter_ctr = 0;  ctr = 1;
    
    record_U_g= 1;
    record_U_a= 1;
    
    record_r22g = 1;
    record_r11g = 1;
    record_n21g = 1;
    
    record_r22a = 1;
    record_r11a = 1;
    record_n21a = 1;
    
end

%%%% specify some of the mainloop control parameters %%%%
iter_per_rt = round(dat.T_R/dat.dt);
checkptIter = iter_per_rt*100; %make a checkpoint every 100RTs.
% sim_settings.simRT = 200;
tEnd = sim_settings.simRT*dat.T_R; % end time in tps
plotCtr = sim_settings.plotCtr; %% set on how many iterations should the program plot the intensity

%simulation info storage arrays -> preallocate
recordingduration = sim_settings.recordRT*dat.T_R; % how many ps should we record the pulse for

iterperrecord = 1;
recordingiter  = round(recordingduration/iterperrecord/dat.dt);
padsize = double(recordingiter-length(record_U_g));

%preallocate memory for optical field, curr density voltage wave, etc.
%strorage...
record_U_g= padarray(record_U_g,padsize,'post');
record_U_a = padarray(record_U_a,padsize,'post');

%store population info
record_r22g = padarray(record_r22g,padsize,'post');
record_r11g = padarray(record_r11g,padsize,'post');
record_r22a = padarray(record_r22a,padsize,'post');
record_r11a = padarray(record_r11a,padsize,'post');

record_n21g = padarray(record_n21g,padsize,'post');
record_n21a = padarray(record_n21a,padsize,'post');

if(strcmp(dat.dtype,'single'))
    
    record_U_g= single(record_U_g);
    record_U_a = single(record_U_a);
    
    record_n21g = single(record_n21g);
    record_n21a = single(record_n21a);
    
    %store population info
    record_r22g = single(record_r22g);
    record_r11g = single(record_r11g);
    record_r22a = single(record_r22a);
    record_r11a = single(record_r11a);
    
end

info.settings = sim_settings;
info.cavity = 'RING-SITMODELOCKING';
info.Ltot = dat.Ltot;
info.N = sim_settings.N;
info.SIMTYPE = ['SIT modelocking. Abs length: ' num2str(sim_settings.La) 'mm; Gain length: ' num2str(sim_settings.Lg) ' mm'];
r22 = zeros(sim_settings.N,1,dat.dtype); 
r11 = zeros(sim_settings.N,1,dat.dtype); 

while(dat.t< tEnd)
    
    
    %%%%% end of setting up the TL PARAMS %%%%%
    
    if(mod(iter_ctr+1,checkptIter) == 0 )
        checkptname = ['CHCKPT_' sim_settings.name '_SIT_' num2str(sim_settings.N) '_FP'];
        save(checkptname);
    end
    
    %%plot some of the results if neeed ariseth :D
    if(mod(iter_ctr,500) == 0)
        clc;
        info.iter_ctr = iter_ctr;
        info.RT = dat.t/dat.T_R;
        intensity = dat.U.*conj(dat.U) ;
        info.maxInt  =  max(intensity);
        printINFO(info);
        
        
        r22(dat.gain_start:dat.gain_end) = dat.r22g;
        r11(dat.gain_start:dat.gain_end) = dat.r11g;
        r22(dat.abs_start:dat.abs_end) = dat.r22a;
        r11(dat.abs_start:dat.abs_end) = dat.r11a;
        ax = plotyy(dat.x*1e3,abs(dat.U).^2, dat.x*1e3,r22-r11);
        set(ax(1).Children,'color','b');
        ylabel(ax(1), 'I(t)');
        set(ax(2).Children,'color','r');
        set(ax(2),'YLim',[-1,1]);
        ylabel(ax(2), 'pop. inversion');
        set(ax,{'ycolor'},{'b';'r'})
        title(info.SIMTYPE);
        
        
        getframe;
    end
    
    %%%% obtain the field, field intensity and the total population at position "idx" ...
    %store fields info
    if dat.gain_end-dat.gain_start >0 
        record_U_g(ctr)= dat.U(idx_g);
        record_r22g(ctr) = dat.r22g(idx);
        record_r11g(ctr) = dat.r11g(idx);
        record_n21g(ctr) = dat.n21g(idx);
    end
    if dat.abs_end-dat.abs_start >0 
        record_U_a(ctr)= dat.U(idx_a);
        record_r22a(ctr) = dat.r22a(idx);
        record_r11a(ctr) = dat.r11a(idx);
        record_n21a(ctr) = dat.n21a(idx);
    end
    
    ctr = ctr+1;
    
    dat = stepBloch(sim_settings,dat);
    dat = stepWave(sim_settings,dat);
    dat = updateBloch(sim_settings,dat);
    
    dat.t = dat.t+dat.dt;
    
    iter_ctr = iter_ctr + 1;
end
savename = [sim_settings.name '_N_' num2str(sim_settings.N) '_RING_Overlap1_' num2str(sim_settings.simRT) ];
save(savename);

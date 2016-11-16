clear; clc;close all;
scenariofile = 'RTSLOWLIGHT_PROTOTYPE5-5e16doping.sim';
simfile = 'RTSLOWLIGHT_PROTOTYPE5-5e16doping.sim';

%parse all input files
sim_settings = parseSimParams(scenariofile);
sim_settings = parseSimData(simfile,sim_settings);


init = 1

if(init > 0)
    
    clc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZE DEVICE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%
    dat.Ltot = 2; % mm
    sim_settings.N = 2000;
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
    
    %%%%% Slow light - gain region definition %%%%
    % number of regions
    NR = 100;
    % total length of slowing regions
    L_SD =  dat.Ltot*0.5;
    dW = L_SD/NR;
    
    % length of free space
    L_GAIN = (dat.Ltot-L_SD);
    %width of a single Free space region
    dF = L_GAIN/(NR+1);
    
    % free-space slow-down free-space slow-down free-space
    core_i = zeros(length(dat.x),NR);
    core = 0*dat.x;
    for i = 1:NR
        core_i(:,i) = (dat.x>= (i*dF+(i-1)*dW)  & dat.x<= (i*dF+i*dW));
        core = core + core_i(:,i);
    end
    % interaction area !
    dat.i_core = (dat.x>(dF) & dat.x<=(NR*dF+NR*dW));
    % slow light and gain cores 
    dat.core_SD = core; dat.core_GAIN = 1-core;
    
    %%% SLOW LIGHT SECTION PARAMS %%%
    
    %indices of the upper and lower laser level
    Spin_ = 3; Ex_ = 2; Grnd_ = 1;
    dat.Spin_ = Spin_; dat.Ex_ = Ex_ ; dat.Grnd_  = Grnd_;
    %eigenenergies of the gain hamiltonian -> transform into ang. freq.
    dat.NLVLS = 3;
    HTB = reshape(sim_settings.HTB',dat.NLVLS,dat.NLVLS);
    
    dat.E_s = HTB(Spin_,Spin_)/dat.hbar;
    dat.E_e = HTB(Ex_,Ex_)/dat.hbar;
    dat.E_g = HTB(Grnd_,Grnd_)/dat.hbar; % in rad/ps
    dat.E_s = dat.E_e;
    
    dat.O_se = HTB(Spin_,Ex_)/dat.hbar; % in rad/ps
    dat.E_se = dat.E_s-dat.E_e; %rad/ps; 1->3 traisition freq
    dat.E_sg = dat.E_s-dat.E_g; %rad/ps; 1->2 transition freq
    dat.E_eg = dat.E_e-dat.E_g; %rad/ps  3->2 transition freq (optical transition)
    %central field frequency.
    dat.E0 = (dat.E_eg+dat.E_sg)/2;
    
    dat.W = sim_settings.Wmtx;
    NLVLS = sqrt(length(dat.W)); % nr of levels to consider !
    %reshape into a matrix
    dat.W = reshape(dat.W,NLVLS,NLVLS).';
    for j = 1:NLVLS
        dat.W(j,j) = 0;
    end
    dat.G = zeros(NLVLS,1);
    zeroFCTR = 1;
    
    %slow down layer lifetimes
    dat.G(Spin_,1) = sum(dat.W(Spin_,:));
    dat.G(Ex_,1) = sum(dat.W(Ex_,:));
    dat.G(Grnd_,1) = sum(dat.W(Grnd_,:));
   
    
 
    gamma_se = (  1/2*(dat.G(Spin_) +  dat.G(Ex_))  ); %% dephsing of the resonant tunneling transition
    gamma_eg = (  1/2*(dat.G(Ex_) +  dat.G(Grnd_))  ); % dephasing of the optical transision...
    gamma_sg = (  1/2*(dat.G(Spin_)+ dat.G(Grnd_))  ); % dephasing of the latest transition

    dat.dE_se = -1i*dat.E_se - gamma_se; %
    dat.dE_eg = +1i*(dat.E0 - dat.E_eg) - gamma_eg; %
    dat.dE_sg = +1i*(dat.E0 - dat.E_sg)- gamma_sg; %
    
    
    
    
    %slow down layer Carrier density and trace!
    zUL_SD = sim_settings.zUL;
    dat.Ncarriers_SD = sim_settings.dN*(100^3)*sim_settings.Ld/sim_settings.Lp; % cm^-3 --> m^-3; carrier density
    dat.trace_rho_SD = ((dat.E0*1E12*dat.Ncarriers_SD*sim_settings.Overlap* ...
        ((zUL_SD*1E-9*Constants('q0'))^2))/(Constants('eps0')*sim_settings.nTHz*Constants('c')*Constants('hbar')))/(1/(sim_settings.lch*sim_settings.tch));
    %%% END SLOW LIGHT SECTION PARAMS %%%
    
    %%% GAIN SECTION PARAMS %%%%
    %resonance energy -> make it the same as the slow light section
    %% gain recovery and dephasing times
    dat.T2 = 2.35;
    dat.T1 = 10;
    dat.dE21 = -1/dat.T2;
    
    %doping density
    avg_dN_g = 2.3e17; % in cm^-3
    zUL_GAIN = 4.5;
    dat.Ncarriers_GAIN = avg_dN_g*(100^3); % cm^-3 --> m^-3; carrier density
    Overlap = sim_settings.Overlap;
    dat.trace_rho_GAIN = ((dat.E0*1E12*dat.Ncarriers_GAIN*Overlap*...
        ((zUL_GAIN*1E-9*Constants('q0'))^2))/(Constants('eps0')* ...
        sim_settings.nTHz*Constants('c')*Constants('hbar')))...
        /(1/(sim_settings.lch*sim_settings.tch));
    dat.r22_0 = 1;
    dat.r11_0 = 0;
    
    %%% END GAIN SECTION PARAMS %%%
    
    %the varying dipole ratio.
    dat.dipR = ones(sim_settings.N,1,dat.dtype)*zUL_GAIN/zUL_SD.*dat.core_GAIN ...
        + dat.core_SD;
    dat.dipR_inv = 1./dat.dipR;
    
    dat = makeMaxwellVars(sim_settings,dat);
    dat = makeBlochVars(sim_settings,dat);
    
    dat.t = dat.dt;
    % the index of the predefined point iun the
    % gain and absorption media where to sample the feild at!
    
    idx = 1;
    iter_ctr = 0;  ctr = 1;
    
    U_in= 1;
    U_out= 1;
    
end

%%%% specify some of the mainloop control parameters %%%%
iter_per_rt = round(dat.T_R/dat.dt);
checkptIter = iter_per_rt*100; %make a checkpoint every 100RTs.
sim_settings.simRT = 90;
tEnd = sim_settings.simRT*dat.T_R; % end time in tps
plotCtr = sim_settings.plotCtr; %% set on how many iterations should the program plot the intensity
interpCtr = 100; %set how often to interpolate the energies, scattering rates, dipole elements etc.

%simulation info storage arrays -> preallocate
recordingduration = sim_settings.recordRT*dat.T_R; % how many ps should we record the pulse for

iterperrecord = 1;
recordingiter  = round(recordingduration/iterperrecord/dat.dt);
padsize = double(recordingiter-length(U_in));

%preallocate memory for optical field, curr density voltage wave, etc.
%strorage...
U_in= padarray(U_in,padsize,'post');
U_out = padarray(U_out,padsize,'post');

%store population info

if(strcmp(dat.dtype,'single'))
    
    U_in= single(U_in);
    U_out = single(U_out)*(dat.t>100);
    
end

info.settings = sim_settings;
info.cavity = 'SLOWLIGHT-ABSORPTION-GAIN';
info.Ltot = dat.Ltot;
info.N = sim_settings.N;
info.SIMTYPE = ['SLOW LIGHT WITH GAIN SECTION. Slow light region length: ' ...
    num2str(L_SD) 'mm; Gain region length: ' num2str(L_GAIN) ' mm'];
r22 = zeros(sim_settings.N,1,dat.dtype);
r11 = zeros(sim_settings.N,1,dat.dtype);
vid_ctr = 1;
while(dat.t< tEnd)
    
    
    %%%%% end of setting up the TL PARAMS %%%%%
    
    if(mod(iter_ctr+1,checkptIter) == 0 )
        checkptname = ['CHCKPT_' sim_settings.name '_SIT_' num2str(sim_settings.N) '_FP'];
        save(checkptname);
    end
    
    %%plot some of the results if neeed ariseth :D
    if(mod(iter_ctr,500) == 0)
        
        clc;
        gain = ((dat.T2*1e-12*dat.E0*1E12*dat.Ncarriers_GAIN*Overlap*...
        ((zUL_GAIN*1E-9*Constants('q0'))^2))/(Constants('eps0')*sim_settings.nTHz*Constants('c')*Constants('hbar')));
        display(['gain = ' num2str(gain) ])
    
        info.iter_ctr = iter_ctr;
        info.RT = dat.t/dat.T_R;
        intensity = dat.U.*conj(dat.U) ;
        info.maxInt  =  max(intensity);
        printINFO(info);
        
        
        %%%%%%%%%%%%%%%%%
%         subplot(2,1,1)
        plot(dat.x,abs(dat.U).^2,'b','Linewidth',2.0);
        xlabel('x (mm)');ylabel('Intensity (a.u.)');
        hold on;
        fobj = fill(dat.x,dat.core_SD*dat.ampl^2,'g'); 
        set(fobj,'EdgeColor','none','FaceAlpha',0.3);
        hold off;
%         title(info.SIMTYPE);
        title(['t = ' num2str(dat.t) ' (ps)'])
%         subplot(2,1,2);
%         r22 = dat.r_ee.*dat.core_SD+dat.r22.*dat.core_GAIN;
%         r11 = dat.r_gg.*dat.core_SD+dat.r11.*dat.core_GAIN;
%         plot(dat.x,r22-r11,'-r');
%         legend('\Delta(x,t)');
        %%%%%%%%%%%%%%%%%
        M_ovie(vid_ctr) = getframe;
%         vid_ctr = vid_ctr+1;
    end
    %%%% obtain the field, field intensity and the total population at position "idx" ...
    
    %store fields info
    U_in(ctr)= dat.U(1);
    U_out(ctr)= dat.U(end);
    ctr = ctr+1;

    dat = stepBloch(sim_settings,dat);
    dat = stepWave(sim_settings,dat);
    dat = updateBloch(sim_settings,dat);
    
    
    dat.t = dat.t+dat.dt;
    iter_ctr = iter_ctr + 1;
    
end
% v = VideoWriter('message_transmissionx2.avi');
% open(v)
% for frame_num = 1:length(M_ovie)
% writeVideo(v,M_ovie(frame_num));
% end
% close(v);
% savename = [sim_settings.name '_N_' num2str(sim_settings.N) '_' num2str(sim_settings.simRT) ];
% save(savename);

%% clean clear and close stuff
close all; clc; 

convention  = 'physics';

%determine the transform function according to the assumed convention
if(strcmp(convention,'physics'))
    transform = @ifft;
    itransform = @fft; 
    cs = 1;
else
    transform = @fft;
    itransform = @ifft;
    cs = -1;
end

rt_start = 500; rt_end =995;
Dt = dt*iterperrecord; 

iter_per_rt = round(T_R/Dt);
envelope = E_p(rt_start*iter_per_rt+1:rt_end*iter_per_rt)+E_m(rt_start*iter_per_rt+1:rt_end*iter_per_rt);
% envelope = E_p(rt_start*iter_per_rt+1:rt_end*iter_per_rt);

params.N = length(envelope); 
params.dt = Dt; 
params.iterperrec = iterperrecord; 
params.dv = f_R;
params.E0 = E0; 
params.convention = convention;
params.cs = cs; 
params.Ntau = 100*iter_per_rt; 

dat = calcswiftNEW(envelope,params); 

limI = dat.lags > -10*T_R & dat.lags < 10*T_R;
dat.SWIFTS = dat.SWIFTS(limI,:);  tms = dat.lags(limI); 

THz_range   = [1,7];
ZERO_PAD_DF = 100e-6;          % in THz


% x = load(data);
% iqn = [getfield(x,['v',num2str(iqn_channels(1))]).',...
%        getfield(x,['v',num2str(iqn_channels(2))]).',...
%        getfield(x,['v',num2str(iqn_channels(3))]).'];
%    
iqn = dat.SWIFTS;
fs = 1/Dt;
df   = fs/size(iqn,1);
% Dt   = 1/(size(iqn,1)*df);
Niqn    = size(iqn,1);

iqn = iqn./repmat(polyval(polyfit([1:Niqn]',iqn(:,3),1),[1:Niqn]'),1,3);

A = [[1:Niqn]',ones(Niqn,1)];
iqn = iqn - A*((A'*A)\(A'*iqn)); % demean

%% Determine rep rate and IQ calibration factor
N_zp     = ceil(1/(ZERO_PAD_DF*Dt));
fs_zp    = 1/(N_zp*Dt)*[0:N_zp-1];
Fs_zp    = transform(iqn.*repmat(hanning(size(iqn,1)),1,3),N_zp);
[~,ml]= max(abs(Fs_zp(fs_zp<THz_range(1),1)));
rep_rate = 1000*params.dv         % in GHz

gis = find(fs_zp<THz_range(1));
C = (Fs_zp(gis,1)'*Fs_zp(gis,1))\(Fs_zp(gis,1)'*Fs_zp(gis,2))

iqnpm = [iqn, (C*iqn(:,1) - iqn(:,2))/(1i*imag(C)), ...
       (conj(C)*iqn(:,1) - iqn(:,2))/(-1i*imag(C))];

%  iqnpm(:,4) = conj(iqnpm(:,4));
%  iqnpm(:,5) = conj(iqnpm(:,5));

   
iqnpm = [iqnpm,abs(iqnpm(:,4))];

%% Find peaks of normal spectrum. Deconvolve data.
Fs_zp    = transform(iqnpm(:,3).*hanning(Niqn),N_zp);
gis = find(fs_zp>THz_range(1) & fs_zp<THz_range(2));
[pks,pk_fs] = findpeaks(abs(Fs_zp(gis)),fs_zp(gis)','MinPeakDistance',rep_rate/1000*.75);

ts    = Dt *[0:size(iqnpm,1)-1]';
exps  = exp(-1i*2*pi*repmat(pk_fs.',length(ts),1).*repmat(ts,1,length(pk_fs)));
we = repmat(sqrt(hanning(Niqn)),1,size(exps,2)); wp = repmat(sqrt(hanning(Niqn)),1,size(iqnpm,2));
pk_cs = ((exps.*we)'*(exps.*we))\((exps.*we)'*(iqnpm.*wp));


%%
Spp = sqrt(abs(pk_cs(1:end-1,3).*pk_cs(2:end,3)));
Sp = pk_cs(1:end-1,4);
Spp = Spp * ((Spp'*Spp)\(Spp'*abs(Sp))); %abs(Sp(ml))/Spp(ml);

%%
fs    = 1/(Niqn*Dt)*[0:Niqn-1];
dfigure;
set(gcf,'Position',[559 133 765 919]); movegui(gcf,'center');
subplot(4,1,1); plot(ts,iqn);
title('De-meaned IQN Data'); xlabel('Time (ps)'); ylabel('DAQ Signal');
legend('I','Q','N');
subplot(4,1,2); semilogy(fs,abs(transform(iqnpm(:,[1,2,3]).*repmat(hanning(Niqn),1,3)))); xlim(THz_range);
title('FT of IQN (Hanning-windowed)'); xlabel('Frequency (THz)'); ylabel('Intensity (a.u.)');
legend('I','Q','N');
subplot(4,1,3); semilogy(fs,abs(transform(iqnpm(:,[4,5,3]).*repmat(hanning(Niqn),1,3)))); xlim(THz_range);
title('FT of PMN (Hanning-windowed)'); xlabel('Frequency (THz)'); ylabel('Intensity (a.u.)');
legend('P','M','N');
subplot(4,1,4); semilogy(pk_fs(1:end-1),abs(Sp),pk_fs(1:end-1),Spp);
title('Deconvolved coefficients'); xlabel('Frequency (THz)'); ylabel('Intensity (a.u.)');
legend('c_+^{(i)}','(c_0^{(i)}c_0^{(i+1)})^{1/2}');
hold all; plot(pk_fs,abs(pk_cs(:,end)));

%%
output.ts = ts;
output.fs = fs';
output.iqnpm = iqnpm;
output.C = C;
output.pk_fs = pk_fs;
output.pk_cs = pk_cs;
% 
% x = SWIFTS_TD_alt(output,[2,3],[],[])
% xlim([2.5,4.5]); ylim([0,1.2]);
% xlabel('Frequency (THz)'); ylabel('g(\omega)'); title('Degree of coherence');


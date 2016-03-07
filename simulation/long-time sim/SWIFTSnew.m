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

apodizefunc = @hanning; 


rt_start = 500; rt_end =995;
Dt = dt*iterperrecord; 

iter_per_rt = round(T_R/Dt);
envelope = E_p(rt_start*iter_per_rt+1:rt_end*iter_per_rt)+E_m(rt_start*iter_per_rt+1:rt_end*iter_per_rt);
% envelope = E_p(rt_start*iter_per_rt+1:rt_end*iter_per_rt);
envelope = reshape(envelope,length(envelope),1);

Npts = length(envelope);
NFFT = Npts;
win = apodizefunc(Npts);

P_t = abs(envelope).^2; Y_note =transform(P_t.*win,NFFT);
% truncate Y_note, because P_t is a real signal.
f_note = 1/Dt*[0:NFFT-1]/NFFT; 
[pks,freqs] = findpeaks(abs(Y_note(f_note > f_R/2)),f_note((f_note > f_R/2)),'MinPeakDistance',f_R*0.75);
dv = freqs(1); 


params.N = length(envelope); 
params.dt = dt; 
params.iterperrec = iterperrecord; 
params.dv = dv;
params.E0 = E0; 
params.convention = convention;
params.cs = cs; 
params.Ntau = 100*iter_per_rt; 

dat = calcswift_II(envelope,params); 


%downsample the swifts 
fs_dwn = 20; %in THz 
fs = 1/Dt; 
sample_freq = round(fs/fs_dwn); 
%take only 10 RTs and downsample at the same time! 

%downsample swifts and take data only for delay = +-5 RT. 
SWIFTS = dat.SWIFTS(1:sample_freq:end,:);
tms = dat.lags(1:sample_freq:end).';
limI = tms > -5*T_R & tms < 5*T_R;
SWIFTS = SWIFTS(limI,:);  tms = tms(limI); 

% restrict SWIFTS to frequency of interest only: 

%at this point we have the downsampled coeffs, frequencies and SWIFT
%interferogram data! ToDo : find the peaks, least square fit, and compute
%the spectrum.


coeffs = transform(SWIFTS.*repmat(apodizefunc(size(SWIFTS,1)),1,size(SWIFTS,2)),NFFT); 
f = fs/sample_freq/size(coeffs,1)*[0:NFFT-1].';
limI = f > 1 & f < 7;  coeffs = coeffs(limI,:);  f = f(limI,:); df = f(2)-f(1); 

[coeffs_,locs]  = findpeaks(abs(coeffs(:,1)),f,'MinPeakDistance',0.75*dv); 
Fmtx= repmat(locs.',size(tms,1),1); Tmtx = repmat(tms,1,size(locs,1));

Nsig = size(SWIFTS,1); 
we = repmat(sqrt(apodizefunc(Nsig)),1,size(locs,1));
wp = repmat(sqrt(apodizefunc(Nsig)),1,size(SWIFTS,2));
X = exp(-1i*2*pi*Fmtx.*Tmtx);
% C = ((X.*we)'*(X.*we))\((X.*we)'*(SWIFTS.*wp));
C = ((X)'*(X))\((X)'*(SWIFTS));

SWIFTS_deconv = 2*real(X*C);     
tms_deconv = linspace(tms(1),tms(end),size(SWIFTS_deconv,1));

c0 = C(:,1); cp = C(:,3)-1i*C(:,4); 

dfigure;  subplot(1,2,1);
semilogy(locs(1:end-1),sqrt(abs(c0(1:end-1).*c0(2:end))),locs(1:end-1),abs(cp(1:end-1))); 
subplot(1,2,2); 
plot(locs(1:end-1),abs(cp(1:end-1))./sqrt(abs(c0(1:end-1).*c0(2:end))));




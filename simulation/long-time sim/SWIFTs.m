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

%windowing function and length of the fft !
winfunc = @hanning; 
NFFTfunc = @(num) 2^nextpow2(num); 


load('disp_compensated');
c = 0.299792458/3.6; % central freq phase vel.  
T_R = 2*5/c; %round trip time in a 5mm cavity
f_R = 1/T_R; 
E0 = 3.88*2*pi; % central frequency (i.e. energy difference between ULL and LLL); 

Npts = length(envelope); 
P_t = abs(envelope).^2; Y_note =transform(P_t);
% truncate Y_note, because P_t is a real signal.
f_note = 1/Dt*[0:Npts-1]/Npts; 
display('Extracting the RT freq. ...'); 
[pks,freqs] = findpeaks(abs(Y_note(f_note > f_R/2)),f_note((f_note > f_R/2)),'MinPeakDistance',f_R*0.75);
dv = freqs(1); 
%%
params.N = length(envelope); 
params.Dt = Dt; 
params.dv = dv;
params.E0 = E0; 
params.convention = convention;
params.winfunc = winfunc;
params.cs = cs; 

display('Calculating SWIFTs (time domain data) ...');
dat = calcswift_III(envelope,params); 

NFFT = NFFTfunc(1E5);

%downsample the swifts 
fs_dwn = 20; %in THz 
fs = 1/Dt;  sample_freq = round(fs/fs_dwn); 
%take only 10 RTs and downsample at the same time! 

%downsample swifts and take data only for delay = +-5 RT. 
display('Downsampling SWIFTs ...');

SWIFTS = dat.SWIFTS(1:sample_freq:end,:);
tms = dat.lags(1:sample_freq:end).';
limI = tms > -5*T_R & tms < 5*T_R;
SWIFTS = SWIFTS(limI,:);  tms = tms(limI); 

%at this point we have the downsampled coeffs, frequencies and SWIFT
%interferogram data! ToDo : find the peaks, least square fit, and compute
%the spectrum.


coeffs = transform(SWIFTS.*repmat(winfunc(size(SWIFTS,1)),1,size(SWIFTS,2)),NFFT); 
f = fs/sample_freq/size(coeffs,1)*[0:NFFT-1].';
limI = f > 1 & f < 7;  coeffs = coeffs(limI,:);  f = f(limI,:); df = f(2)-f(1); 

display('Identifying the spectral freq. peaks ...');
[coeffs_,locs]  = findpeaks(abs(coeffs(:,1)),f,'MinPeakDistance',0.75*dv); 
Fmtx= repmat(locs.',size(tms,1),1); Tmtx = repmat(tms,1,size(locs,1));

display('Least square''s fitting the  coeffs. ...');
Nsig = size(SWIFTS,1); 
% do not really udnerstand what are these we and wp functions for?
we = repmat(sqrt(winfunc(Nsig)),1,size(locs,1));
wp = repmat(sqrt(winfunc(Nsig)),1,size(SWIFTS,2));

X = exp(-1i*2*pi*Fmtx.*Tmtx);
C = ((X)'*(X))\((X)'*(SWIFTS));

% C = ((X.*we)'*(X.*we))\((X.*we)'*(SWIFTS.*wp));

SWIFTS_deconv = 2*real(X*C);     
tms_deconv = linspace(tms(1),tms(end),size(SWIFTS_deconv,1));

c0 = C(:,1); cp = C(:,3)-1i*C(:,4); 

figure;  
subplot(1,2,1);
semilogy(locs(1:end-1),sqrt(abs(c0(1:end-1).*c0(2:end))),locs(1:end-1),abs(cp(1:end-1))); 
xlabel('Freq. (THz)'); legend('spectrum product','correlation');

subplot(1,2,2); 
plot(locs(1:end-1),abs(cp(1:end-1))./sqrt(abs(c0(1:end-1).*c0(2:end))));
xlabel('Freq. (THz)'); ylabel('coherence');



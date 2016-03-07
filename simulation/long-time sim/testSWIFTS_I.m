
transform = @ifft;
itransform = @fft; 
    
iqn_channels = [4,5,6];
x = load('Data at 0.9002 A, 50 K.mat');
SWIFTS = [getfield(x,['v',num2str(iqn_channels(1))]).',...
       getfield(x,['v',num2str(iqn_channels(2))]).',...
       getfield(x,['v',num2str(iqn_channels(3))]).'];

   
%downsample the swifts 
fs = x.d.f_THz*size(SWIFTS,1); 
df = x.d.f_THz; 
f_R = 6.799992514731420/1000; 
T_R = 1/(f_R);
Dt   = 1/(size(SWIFTS,1)*df);
lags = Dt*size(SWIFTS,1)*linspace(-1/2,1/2,size(SWIFTS,1));

fs_dwn = 20; %in THz 
sample_freq = round(fs/fs_dwn); 
%take only 10 RTs and downsample at the same time! 

%downsample swifts and take data only for delay = +-5 RT. 
SWIFTS = SWIFTS(1:sample_freq:end,:);
tms = lags(1:sample_freq:end).';
limI = tms > -5*T_R & tms < 5*T_R;
SWIFTS = SWIFTS(limI,:);  tms = tms(limI); 


% restrict SWIFTS to frequency of interest only: 

%at this point we have the downsampled coeffs, frequencies and SWIFT
%interferogram data! ToDo : find the peaks, least square fit, and compute
%the spectrum.

NFFT = 1E6;

coeffs = transform(SWIFTS.*repmat(hanning(size(SWIFTS,1)),1,size(SWIFTS,2)),NFFT); 
f = fs/sample_freq/size(coeffs,1)*[0:NFFT-1].';

gis = find(f<1);
SWIFTS = [SWIFTS, SWIFTS(:,1) +1i*SWIFTS(:,2),SWIFTS(:,1) - 1i*SWIFTS(:,2)];


limI = f > 1 & f < 7;  coeffs = coeffs(limI,:);  f = f(limI,:); df = f(2)-f(1); 

[coeffs_,locs]  = findpeaks(abs(coeffs(:,3)),f,'MinPeakDistance',0.75*f_R); 
Fmtx= repmat(locs.',size(tms,1),1); Tmtx = repmat(tms,1,size(locs,1));

Nsig = size(SWIFTS,1); 
we = repmat(sqrt(hanning(Nsig)),1,size(locs,1));
wp = repmat(sqrt(hanning(Nsig)),1,size(SWIFTS,2));
X = exp(-1i*2*pi*Fmtx.*Tmtx);
C = ((X.*we)'*(X.*we))\((X.*we)'*(SWIFTS.*wp));

SWIFTS_deconv = 2*real(X*C);     
tms_deconv = linspace(tms(1),tms(end),size(SWIFTS_deconv,1));

c0 = C(:,3); cp = C(:,4); 

subplot(1,2,1);
plot(locs(1:end-1),sqrt(abs(c0(1:end-1).*c0(2:end))),locs(1:end-1),abs(cp(1:end-1))); 
subplot(1,2,2); 
plot(locs(1:end-1),abs(cp(1:end-1))./sqrt(abs(c0(1:end-1).*c0(2:end))));




function data = calcswiftdeconvolution( env, params )


N = params.N; Dt = params.dt*params.iterperrec;  E0 = params.E0; Ntau = params.Ntau;  dv = params.dv;
env = reshape(env,N,1);
convention = params.convention;

%determine the transform function according to the assumed convention
if(strcmp(convention,'physics'))
    transform = @ifft;
    cs = 1;
else
    transform = @fft;
    cs = -1;
end

apodizefunc = @hanning;
NFFTfunc = @(elem) 2^nextpow2(elem);
% apodizefunc = @(elem)ones(elem,1);
% NFFTfunc = @(elem) elem;


T_tot = N*Dt; tms = T_tot*linspace(0,1,N)';fs = 1/Dt;

%E_t and E_t-dv
E_t = real(env.*exp(-cs*1i*E0*tms));
E_t_shift = real(env.*exp(-cs*1i*E0*tms).*exp(cs*1i*2*pi*dv*tms));

%take the normal interferogram with a window of param.Ntau*Dt;
[S0,lags0] = xcorr(E_t,E_t,Ntau);
S0 = S0.*hanning(length(S0));
% analytic signal
[SP,lagsP] = xcorr(E_t.*exp(1i*2*pi*dv*tms),E_t,Ntau);
SP = SP.*hanning(length(SP));

f = fs*linspace(0,1/2,length(S0));

BW = 2; 
lim_idx = (f>=E0/2/pi-BW/2) & (f <=E0/2/pi +BW/2);
f_ = f(lim_idx); 

lags_= linspace(-1/2,1/2,length(f_))*Dt; 


X = zeros(length(lags_));
for p = 1:length(lags_)
    tau = lags_(p);
    X(p,:) = exp(-cs*1i*2*pi*f_*tau);
end

c0 = X'*X\X'*S0;
S0_lsqr = X*c0;
err0 = (S0-S0_lsqr)'*(S0-S0_lsqr);


cp  = X'*X\X'*SP;
SP_lsqr = X*cp;
errP = (SP-SP_lsqr)'*(SP-SP_lsqr);

data.E_t =E_t; data.tms =tms;
data.S0 = S0 ; data.lags0= lags0;  data.SP = SP; data.lagsP = lagsP;
data.lags_ = lags_;
data.S0_lsqr = S0_lsqr; data.SP_lsqr = SP_lsqr;
data.err0 =  err0 ; data.errP = errP;


%freq domain data!
data.f = f; data.f_ = f_; 
data.c0 = c0 ; data.cp = cp;

end


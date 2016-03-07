function data = calcswift_I( envelope, params )


N = params.N; Dt = params.dt*params.iterperrec;  E0 = params.E0; Ntau = params.Ntau;  dv = params.dv;
envelope = reshape(envelope,N,1);
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
NFFTfunc = @(elem) 2^nextpow2(2*elem); 
% apodizefunc = @(elem)ones(elem,1); 
% NFFTfunc = @(elem) elem; 



%iterations per round trip
tms = Dt*[0:N-1]'; fs = 1/Dt;



if strcmp(params.complex,'yes')
    
    E_t = (envelope.*exp(-cs*1i*E0*tms));
    
    % first do the fourier transform of E
    NFFT = NFFTfunc(N); win = apodizefunc(N);
    nconst = normalizefourier(NFFT,win,transform); 
    E_w = transform(E_t.*win,NFFT)/nconst;  f_w = fs*linspace(0,1,length(E_w));
    
    %take the normal interferogram with a window of param.Ntau*Dt;
    [S0,lags] = xcorr((E_t),E_t,Ntau);   
    Ninterf = length(S0);

    %calculate normalization constant
    NFFT = NFFTfunc(Ninterf);
    win = apodizefunc(Ninterf); 
    nconst = normalizefourier(NFFT,win,transform); 
    
    P0 = transform(S0.*win,NFFT)/nconst;   f0 = fs*linspace(0,1,length(P0));
    E_t_shift = E_t.*exp(cs*1i*2*pi*dv*tms); [S1,lags] =  xcorr((E_t_shift),E_t_shift,Ntau);
    P1 = transform(S1.*win,NFFT)/nconst;
    
    %in-phase and in-quadrature components
    [SI,lags] = xcorr(conj(E_t.*cos(2*pi*dv*tms)),E_t,Ntau); [SQ,lags] = xcorr(conj(E_t.*sin(2*pi*dv*tms)),E_t,Ntau);
    
    %results: swift time domain signal
    SP = SI+1i*SQ;
    
    %swift frequency domain signal
    SWIFT = transform(SP.*win,NFFT)/nconst;
    
else
    
    %E_t and E_t/dv 
    E_t = real(envelope.*exp(-cs*1i*E0*tms));
    E_t_shift = real(envelope.*exp(-cs*1i*E0*tms).*exp(cs*1i*2*pi*dv*tms));
    
    % E_w
    NFFT = NFFTfunc(N); win = apodizefunc(N);
    nconst = normalizefourier(NFFT,win,transform); 
    E_w =transform(E_t.*win,NFFT)/(nconst/2); E_w = E_w(1:NFFT/2+1); f_w = fs*linspace(0,1/2,length(E_w));
    
    %take the normal interferogram with a window of param.Ntau*Dt;
%     [S0,lags0] = xcorr(E_t,E_t,Ntau);
%     Ninterf = length(S0); 
%     NFFT =NFFTfunc(Ninterf); 
%     win = apodizefunc(Ninterf); 
%     nconst = normalizefourier(NFFT,win,transform); 
%     P0 = transform(S0.*win,NFFT)/(nconst/2);  P0 = P0(1:length(P0)/2+1); f0 = fs*linspace(0,1/2,length(P0));
%     
%     [S1,lags1] =  xcorr(E_t_shift,E_t_shift,Ntau);
%     P1 = transform(S1.*win,NFFT)/(nconst/2);  P1 = P1(1:length(P1)/2+1);
%     
%     %in-phase and in-quadrature components
%     [SI,lagsI] = xcorr(E_t.*cos(2*pi*dv*tms),E_t,Ntau); [SQ,lagsQ] = xcorr(E_t.*sin(2*pi*dv*tms),E_t,Ntau);
%     
%     %results: swift time domain signal
%     SP = SI+cs*1i*SQ;
%     %swift frequency domain signal
%     SWIFT = transform(SP.*win,NFFT)/(nconst/2); SWIFT = SWIFT(1:length(SWIFT)/2+1);


    [S0,lags] = xcorr(E_t,E_t,Ntau);
    S0 = flipud(S0);
    
    Ninterf = length(S0); 
    NFFT =NFFTfunc(Ninterf); 
    win = apodizefunc(Ninterf); 
    nconst = 1;%normalizefourier(NFFT,win,transform); 
    P0 = transform(S0.*win,NFFT)/(nconst); f0 = fs*[0:NFFT-1]/NFFT;
    
    [S1,lags] =  xcorr(E_t_shift,E_t_shift,Ntau);
    S1flip = flipud(S1);
    
    P1 = transform(S1.*win,NFFT)/(nconst);  
    
    %in-phase and in-quadrature components
    [SI,lags] = xcorr(E_t.*cos(2*pi*dv*tms),E_t,Ntau); 
    SI = flipud(SI);
    [SQ,lags] = xcorr(E_t.*sin(2*pi*dv*tms),E_t,Ntau);
    SQ = flipud(SQ);
    

    %results: swift time domain signal
    SP = SI-cs*1i*SQ;
    %swift frequency domain signal
    SWIFT = transform(SP.*win,NFFT)/(nconst); 

end

%time domain data:
data.E_t = E_t; data.times = tms;
data.S0 = S0; data.lags = lags*Dt;  data.S1 = S1; 
data.SI = SI;  data.SQ = SQ;

%freq domain data
data.E_w = E_w; data.f_w = f_w; data.P0 = P0; data.f0 = f0;
data.P1 = P1; data.SWIFT = SWIFT;

end


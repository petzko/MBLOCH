function [ Interf,lags] = InterferogramFast( f,g,dt,varargin)
    
    


    skipctr = 1; options = {'no'};
%     if (length(varargin) >= 1)
%         skipctr = varargin{1};
%     end
%     if (length(varargin) >= 2)
%         options = varargin{2};
%     end
%     
   
   

    N = length(f);
    f_1 = reshape(f,N,1); 
    g_1 = reshape(g,N,1); 
    
    NFFT = N;

    if (strcmp('hamming',options{1}))
        f_1 = f_1.*hamming(N); % window the signal!
        g_1 = g_1.*hamming(N);
        NFFT = 2^nextpow2(N);
    elseif (strcmp('chebwin',options{1}))
        f_1 = f_1.*hamming(N); % window the signal!
        g_1 = g_1.*hamming(N);
        NFFT = 2^nextpow2(N);
    elseif (strcmp('barthannwin',options{1}))
       f_1 = f_1.*hamming(N); % window the signal!
        g_1 = g_1.*hamming(N);
        NFFT = 2^nextpow2(N);
    elseif (strcmp('flattopwin',options{1}))
        f_1 = f_1.*hamming(N); % window the signal!
        g_1 = g_1.*hamming(N);
        NFFT = 2^nextpow2(N);
    end

    T = NFFT*dt*skipctr;
    Interf = ifft(conj(fft(f_1,NFFT)).*fft(g_1,NFFT));
    Interf = fftshift(Interf);
    lags = T*(linspace(-1/2,1/2,NFFT));

end
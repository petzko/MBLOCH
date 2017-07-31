function [ X,f ] = mydft2(x,dt)
     % Assumes x is an N x 1 column vector
     % and Fs is the sampling rate.
     dt = dt*1e-12;
     N = length(x);
     
     t = dt*(0:N-1)';
     
     NFFT = 2^nextpow2(N);
     win = hanning(N);
     win = reshape(win,size(x));
     % find the normalization constant
     Ynrom = fft(ones(size(x)).*win,NFFT);  
     nconst = abs(Ynrom(1));
     if isreal(x)
         nconst = nconst/2; 
     end
     X = fft(x.*win,NFFT)/nconst;
%      X = X(1:NFFT/2);
%      f = 1/dt*[0:NFFT/2-1]/NFFT;
    f = 1/dt*[0:NFFT/2-1,-NFFT/2:-1]/NFFT;
 
%      figure;
%      plot(f,normc(abs(X)));
end
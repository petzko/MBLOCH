function [ X,f,t ] = mydft(x,dt)
     % Assumes x is an N x 1 column vector
     % and Fs is the sampling rate.
     dt = dt*1e-12;
     N = length(x);
     
     t = dt*(0:N-1)';
     
     NFFT = 2^nextpow2(N);
     win = hanning(N);
     win = reshape(win,size(x));
     
     X = fft(x.*win,NFFT)/NFFT;
     X = X(1:NFFT/2);
     f = 1/dt*[0:NFFT/2-1]/NFFT;
  
     figure;
     plot(f,abs(X));
end
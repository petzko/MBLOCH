function [ X,f,t ] = mydft(x,dt)
     % Assumes x is an N x 1 column vector
     % and Fs is the sampling rate.
     dt = dt*1e-12;
     N = size(x,1);
     Fs = 1/dt;
     t = dt*(0:N-1)';
     dF = Fs/N;
     f = dF*(0:N/2-1)';
     X = fft(x)/N;
     X = X(1:N/2);
     X(2:end) = 2*X(2:end);
     figure;
     plot(f,abs(X));
end
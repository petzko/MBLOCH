function [ conv] = convolve( f,g)

    N = length(f);
    NFFT = N;
    conv = ifft(conj(fft(f,NFFT)).*fft(g,NFFT));
    conv = fftshift(conv);

end
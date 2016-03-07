rt_start = 10; rt_end = 20;

Dt = dt*iterperrecord;
iter_per_rt = round(T_R/Dt);
A_t = E_p(rt_start*iter_per_rt+1:rt_end*iter_per_rt);
A_t_orig = E_p_orig(rt_start*iter_per_rt+1:rt_end*iter_per_rt);

A_t = reshape(A_t,[length(A_t),1]);  A_t_orig = reshape(A_t_orig,[length(A_t_orig),1]);  

Npts = length(A_t); tms = Dt*[0:Npts-1].';
 
transform = @ifft;
apodizefunc = @hanning; 
NFFTfunc = @(elem) 2^nextpow2(2*elem); 
NFFT = NFFTfunc(Npts);
win = apodizefunc(Npts); 

E_t = real(A_t.*exp(-1i*E0*tms));
E_t_orig = real(A_t_orig.*exp(-1i*E0*tms));


Y1 = transform(E_t.*win,NFFT); 
Y2 = transform(E_t_orig.*win,NFFT); 

f = 1/Dt/NFFT*[0:NFFT/2-1,-NFFT/2:-1].';
dfigure; 
subplot(2,2,[1 2]); 
plot(f,abs(Y1),f,abs(Y2));xlim([3.2 ,4.5]);
subplot(2,2,3);
plot(tms,E_t);
subplot(2,2,4);
plot(tms,E_t_orig);


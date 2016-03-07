%mode order
df = 0.3; 

f_0 = 10; 

tend = 10; 
dt  = 0.01; 
t = -10
times = 0;
ctr = 0;
S = 0;
phase = 0;

while t < tend
    ctr = ctr +1 ; 
  
    CF = exp(2*pi*1i*f_0*t); 
    
    S(ctr) = 0;
    for m = 1:20
        S(ctr) = S(ctr)+CF*exp(2*pi*1i*m*df*t + 2*pi*rand(1,1))
    end
    
    times(ctr) = t;
    t = t+dt;
    plot(times,abs(S).^2)
    getframe;
    
end


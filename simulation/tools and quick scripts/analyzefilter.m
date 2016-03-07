    Ltot = 5
    x = linspace(0,Ltot,N);
    close all; 
    x_0 = +0.5;tp =0.2;  
    aE_in = @(z,time) exp(-(time-(z-x_0)/c).^2/tp^2);
    U = aE_in(x,0);
    U = wgn(N,1,1); 
    U = U/max(abs(U)); 
    PEAK = max(abs(U));    
    
    %sampling rate; 
    ps = 1/dx; 
    %number of samples 
    Ns = length(x); 
    samples = (1:1:Ns)/Ns; 
  
    df = fir1(80,[0.001,0.1]);
    Un = filter(df,1,U); 
    plot(x,U,x,Un); 
    
    U_k = fft(U); Un_k = fft(Un); 
    figure; 
    plot(samples,abs(U_k),samples,abs(Un_k))
function [ SP,freq,Y01,Y02 ] = SpectrumProduct( E_t,times,dt,skipctr,dv,tch,varargin)

    %%% asuming a real signal.. 
    
    if length(varargin)== 1
        options = varargin{1};
    else
        options = {'no','fftshift'};
    end
    
    
    
    
     Eshift = shiftRealSig(E_t,dt,skipctr,dv,times);

    [I01,lags] = InterferogramFast(E_t,E_t,dt,skipctr,{'no','truncate'}); % to obtain E_w;
    [freq01,Y01] = FourierTransf(I01,dt,0,skipctr,tch,options); % to obtain E_w;
    
    [I02,lags] = InterferogramFast(Eshift,Eshift,dt,skipctr,options); % to obtain E_dw;
    [freq02,Y02] = FourierTransf(I02,dt,0,skipctr,tch,options); % to obtain E_w;
    
    SP = sqrt((abs(Y01)).*(abs(Y02)));

    freq = freq01;

end


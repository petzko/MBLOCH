function [SWIFT,S_I,S_Q,lags] = SWIFTinterferogramFast( E_t,times,dt,skipctr,dv,varargin)


    if length(varargin)>= 1
        options = varargin{1};
    else
        options = {'no','truncate'}; 
    end
    
    N = length(E_t);
    E_t = reshape(E_t,[N 1]);
    T =  N*dt*skipctr;
    
    E_I =  E_t.*cos(2*pi*dv*times);
    E_Q =  E_t.*sin(2*pi*dv*times);
    [S_I,lags] = InterferogramFast(E_I,E_t,dt,skipctr,options);
    [S_Q,lags] = InterferogramFast(E_Q,E_t,dt,skipctr,options);
    SWIFT = S_I - 1i*S_Q;

    


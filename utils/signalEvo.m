function [ max_evo, FWHM_evo ] = signalEvo( x, y, threshold)
    
   % calculates the pulse amplitude and full width at half maximum 
   % evolution of a signal containing a train of pulses y(x) .
   % @param x -the signal dependent varaible -> will be used to measure
   % the width of the pulse
   % @param y - the signal wave-form as a function of x, i.e. y = y(x) 
   % @param threshold a small number needed to separate the pulses. I.e.
   % Noise threshold... 
   
   M = length(y);
   Mx = length(x);
   if( M ~= Mx)
       disp('input signal and time/freq domain do not have equal lengths') ; 
       disp('Please check your data and try again'); 
       return;
   end
   
    
    % factor out the array indices where the signal is above 
    % the noise threshold
    ind = y > threshold; 
    % take the pulse delimiters 
    ind = [0 abs(diff(ind))>0]>0;
    m = length(ind); 
    del_idx = 0;
    ctr = 0; 
    for i = 2:m-1
        if (ind(i) >0)
            ctr = ctr +1; 
            del_idx(ctr) = i; 
        end
    end
    
    L = length(del_idx);
    max_evo = zeros(round(L/2),1); 
    FWHM_evo = zeros(round(L/2),1);
    del_idx;
    % now take the evolution of the pulse
    for i = 1:L/2-3
        max_evo(i) = max(y(del_idx(2*i-1):del_idx(2*i)));
        FWHM_evo(i) = fwhm(x(del_idx(2*i-1):del_idx(2*i)),y(del_idx(2*i-1):del_idx(2*i)));
    end


end


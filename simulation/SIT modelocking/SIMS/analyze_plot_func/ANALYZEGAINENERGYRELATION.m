Nx_gain = size(inversion_data_gain_vs_x,1)
Nx_abs = size(inversion_data_abs_vs_x,1); 

plot(inv_gain_x)
PULSE_ENERGY_of_x = zeros(size(pulse_data_vs_x,1),1);




INVERSION_GAIN_of_x = zeros(size(inversion_data_gain_vs_x,1),1);
dinv_gain = gain_model.rho_u_0-gain_model.rho_l_0;

INVERSION_ABS_of_x = zeros(size(inversion_data_abs_vs_x,1),1);





W_sat = 1/gain_model.T1/gain_model.T2; 
for x_idx = 1:size(pulse_data_vs_x,1)
    
    x_idx
    
    envelope_x = pulse_data_vs_x(x_idx,:);
    tmsdata_x = [-length(envelope_x)/2:length(envelope_x)/2-1].'*dt;
    norm_envelope = envelope_x/max(abs(envelope_x));
    
    [pulse_peaks,locs] = findpeaks(abs(norm_envelope),tmsdata_x,'MinPeakHeight',0.5);
    dlocs = locs(ceil(length(locs)/2)) - locs(ceil(length(locs)/2)-1);
    peak_idx = locs(ceil(length(locs)/2))-dlocs/2 < tmsdata_x  & tmsdata_x < locs(ceil(length(locs)/2)) + dlocs/2;
    peak_env  = envelope_x(peak_idx);
    peak_tm = tmsdata_x(peak_idx);
    PULSE_ENERGY_of_x(x_idx) = 1./(peak_tm(end)-peak_tm(1))*trapz(peak_tm,abs(peak_env).^2);

    
    DURATION = tmsdata_x(end)-tmsdata_x(1);
    if x_idx <= Nx_gain
        INVERSION_GAIN_of_x(x_idx) = 1./DURATION*trapz(tmsdata_x,inversion_data_gain_vs_x(x_idx,:)); 
    else
        DURATION = tmsdata_x(end)-tmsdata_x(1);
        INVERSION_ABS_of_x(x_idx-Nx_gain) = 1./DURATION*trapz(tmsdata_x,inversion_data_abs_vs_x(x_idx-Nx_gain,:)); 
    end
    
end
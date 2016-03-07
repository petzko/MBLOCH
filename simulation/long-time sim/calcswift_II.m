function data = calcswift_II( envelope, params )


N = params.N; Dt = params.dt*params.iterperrec;  E0 = params.E0; Ntau = params.Ntau;  dv = params.dv;
envelope = reshape(envelope,N,1);



%iterations per round trip
tms = Dt*[0:N-1]'; 

%E_t and E_t/dv 
E_t = real(envelope.*exp(-1i*E0*tms));
E_t_shift = real(envelope.*exp(-1i*E0*tms).*exp(+1i*2*pi*dv*tms));

% %normal and freq. shifted interferogram! 
% [S0,lags] = xcorr(E_t,E_t,Ntau); S0 = flipud(S0); 
% [S1,lags] =  xcorr(E_t_shift,E_t_shift,Ntau); S1 = flipud(S1);
% %in-phase and in-quadrature components
% [SI,lags] = xcorr(E_t.*cos(2*pi*dv*tms),E_t,Ntau);  SI = flipud(SI);
% [SQ,lags] = xcorr(E_t.*sin(2*pi*dv*tms),E_t,Ntau);  SQ = flipud(SQ);

[S0,lags] = xcorr(E_t,E_t); S0 = flipud(S0); 
[S1,lags] =  xcorr(E_t_shift,E_t_shift); S1 = flipud(S1);
%in-phase and in-quadrature components
[SI,lags] = xcorr(E_t.*cos(2*pi*dv*tms),E_t);  SI = flipud(SI);
[SQ,lags] = xcorr(E_t.*sin(2*pi*dv*tms),E_t);  SQ = flipud(SQ);

%windowed swifts! 
% SWIFTS = [S0,S1,SI,SQ].*repmat(hanning(length(S0)),1,4); 
SWIFTS = [S0,S1,SI,SQ]; 
% SWIFTS = [SI,SQ,S0].*repmat(hanning(length(S0)),1,3); 
data.SWIFTS = SWIFTS;  data.lags = lags*Dt; 


end


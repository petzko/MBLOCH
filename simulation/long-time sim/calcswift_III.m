function data = calcswift_III( envelope, params )


N = params.N; Dt = params.Dt;  E0 = params.E0;dv = params.dv;
envelope = reshape(envelope,N,1); winfunc = params.winfunc;

%times vector
tms = Dt*[0:N-1]'; 

%E_t and E_t-freq shifted! 
E_t = real(envelope.*exp(-1i*E0*tms));
E_t_shift = real(envelope.*exp(-1i*E0*tms).*exp(+1i*2*pi*dv*tms));

% convolutions
[S0,lags] = xcorr(E_t,E_t); S0 = flipud(S0); 
[S1,lags] =  xcorr(E_t_shift,E_t_shift); S1 = flipud(S1);
%in-phase and in-quadrature components
[SI,lags] = xcorr(E_t.*cos(2*pi*dv*tms),E_t);  SI = flipud(SI);
[SQ,lags] = xcorr(E_t.*sin(2*pi*dv*tms),E_t);  SQ = flipud(SQ);

%windowed swifts! 
% SWIFTS = [S0,S1,SI,SQ].*repmat(winfunc(length(S0)),1,4); 

%unwindowed sfiwts!
SWIFTS = [S0,S1,SI,SQ]; 

data.SWIFTS = SWIFTS;  data.lags = lags*Dt; 


end


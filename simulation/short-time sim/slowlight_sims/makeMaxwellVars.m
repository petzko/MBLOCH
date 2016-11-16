function dat = makeMaxwellVars(settings,dat)

FWHM = 30*2/sqrt(2); %psec
t0 = +50;
sig = FWHM/(2*sqrt(2*log(2))); % std. dev
dat.ampl = 0.001;
gauss = @(t,mu,s) dat.ampl*exp(-(t-mu).^2/(2*s^2));

NBits =  8;
NUM2SEND = 125;
code = de2bi(NUM2SEND,NBits);
dat.message = @(t) 0;

for i = 1:NBits
    dat.message = @(t) dat.message(t) + code(i)*gauss(t,t0+i*4*FWHM,sig)
end



dat.U = zeros(settings.N,1);
if (strcmp(dat.dtype,'single'))
    dat.U = single(dat.U);
end
% approximate magnitude
magn = 0;
%electric field vector
dat.U = magn*( (rand(settings.N,1,dat.dtype)-.5)+ ...
    1i*(rand(settings.N,1,dat.dtype) - .5));
%electric field solver
dat.wave_solver = RNFDSolver(settings.N,dat.dx,+1,dat.c, dat.U);

%position dependent linear losses both in the gain and absorption
%sections
% dat.losses = zeros(settings.N,1,dat.dtype);
dat.losses = dat.core_GAIN;
%position dependent normalization constant (trace)!!!
dat.factor =  -1i*dat.c*dat.trace_rho_SD*dat.core_SD.*ones(settings.N,1,dat.dtype)+...
    -1i*dat.c*dat.trace_rho_GAIN*dat.core_GAIN.*ones(settings.N,1,dat.dtype);


end
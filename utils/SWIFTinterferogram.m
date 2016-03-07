function [ SWIFT,tvec,S_I,S_Q ] = SWIFTinterferogram( E_t,dt,skipctr,dv)


N = length(E_t);
E_t = reshape(E_t,[N 1]);
E_t = E_t/max(abs(E_t)); 

%%%% make sure that we padd S so that length(S)-1 divides 4. 
N_s = N + 4-mod(N-1,4);

diff = N_s-N;
E_t_padded = [zeros(floor(diff/2),1); E_t ; zeros(ceil(diff/2),1)];

N_s = length(E_t_padded);
duration = N_s*dt*skipctr; 
tvec = linspace(-duration/2,duration/2,N_s);

N_t = (N_s-1)/2; k = (N_s-1)/4;
tvec = linspace(-duration/4,duration/4,N_t+1);
SWIFT = zeros(N_t+1,1); %

S_I = SWIFT;S_Q = SWIFT;

for p = 0:N_t

    pp = p+1;
    
    S_I(pp) = 0.5*dt*((E_t_padded(k+1)*E_t_padded(k-p+k+1))*cos(k*dt*2*pi*dv));
    S_Q(pp) = 0.5*dt*((E_t_padded(k+1)*E_t_padded(k-p+k+1))*sin(k*dt*2*pi*dv));
    
    for j = (k+1):(3*k-1)
        S_I(pp) = S_I(pp)+dt*((E_t_padded(j+1)*E_t_padded(j-p+k+1))*cos(j*dt*2*pi*dv));
        S_Q(pp) = S_Q(pp)+dt*((E_t_padded(j+1)*E_t_padded(j-p+k+1))*sin(j*dt*2*pi*dv));
    end
    
    S_I(pp) = S_I(pp) + 0.5*dt*((E_t_padded(3*k+1)*E_t_padded(3*k-p+k+1))*cos((3*k+1)*dt*2*pi*dv));
    S_Q(pp) = S_Q(pp) + 0.5*dt*((E_t_padded(3*k+1)*E_t_padded(3*k-p+k+1))*sin((3*k+1)*dt*2*pi*dv));
    
end
SWIFT= S_I - 1i*S_Q;

end
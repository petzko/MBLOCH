function [ I,times_ac,Et_padded,times_s ] = Interferogram( E_t,dt,skipctr )


N = length(E_t);
E_t = reshape(E_t,[N 1]);
% E_t = E_t/max(abs(E_t)); 

%%%% make sure that we padd S so that length(S)-1 divides 4. 
N_s = N + 4-mod(N-1,4);

diff = N_s-N;
Et_padded = [zeros(floor(diff/2),1); E_t ; zeros(ceil(diff/2),1)];

N_s = length(Et_padded);
duration = N_s*dt*skipctr;
times_s = linspace(-duration/2,duration/2,N_s);

N_t = (N_s-1)/2; k = (N_s-1)/4;
times_ac = linspace(-duration/4,duration/4,N_t+1);
I = zeros(N_t+1,1); %


for p = 0:N_t

    pp = p+1;
    I(pp) = 0.5*dt*(((Et_padded(k+1)*Et_padded(k-p+k+1))));
    for j = (k+1):(3*k-1)
        I(pp) = I(pp)+dt*(((Et_padded(j+1)*Et_padded(j-p+k+1))));
    end
    I(pp) = I(pp) + 0.5*dt*(((Et_padded(3*k+1)*Et_padded(3*k-p+k+1))));
    
end



end
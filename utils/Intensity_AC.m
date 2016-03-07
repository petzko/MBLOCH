function [ AC,times_ac,S_padded,times_s ] = Intensity_AC(S,dt)


N = length(S);
S = reshape(S,[N 1]);
S = S/max(abs(S)); 

%%%% make sure that we padd S so that length(S)-1 divides 4. 
N_s = N + 4-mod(N-1,4);

diff = N_s-N;
S_padded = [zeros(floor(diff/2),1); S ; zeros(ceil(diff/2),1)];

N_s = length(S_padded);
len = N_s*dt;
times_s = linspace(-len/2,len/2,N_s);

N_t = (N_s-1)/2; k = (N_s-1)/4;
times_ac = linspace(-len/4,len/4,N_t+1);
AC = zeros(N_t+1,1); %


for p = 0:N_t

    pp = p+1;
    AC(pp) = 0.5*dt*(abs(S_padded(k+1)*S_padded(k-p+k+1)))^2;
    for j = (k+1):(3*k-1)
        AC(pp) = AC(pp)+dt*(abs(S_padded(j+1)*S_padded(j-p+k+1)))^2;
    end
    AC(pp) = AC(pp) + 0.5*dt*(abs(S_padded(3*k+1)*S_padded(3*k-p+k+1)))^2;
    
end

%normalize... 
AC=AC/max(AC);

end


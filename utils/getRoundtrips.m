function [RT] = getRoundtrips(signal,dt,T_R,skipctr,rtrips,offset)


n = length(rtrips);
Dt = dt*skipctr;
iter = floor(T_R/Dt);
RT = zeros(iter,n);
offsetiter = floor(offset/Dt);
for i =1:n
    rTs = rtrips(i); 
    RT(:,i) = signal(-offsetiter+(rTs*iter+1):-offsetiter+((rTs+1)*iter));
end

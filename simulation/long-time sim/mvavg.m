function Y = mvavg(X,M)

Y = 0*X;
for m = floor(M/2)+1:(length(X)-floor(M/2)); 
    Y(m) = sum(X(m-floor(M/2):m+floor(M/2)))/M;
end

end

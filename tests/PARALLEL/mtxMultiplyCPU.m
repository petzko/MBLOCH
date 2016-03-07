function [ C ] = mtxMultiplyCPU( A,B )
    [N,M] = size(A);
    [K,L] = size(B);
    assert(M == K,'Mtx multiply failed! Matrix dimensions do not match');
    C = A*B;
end


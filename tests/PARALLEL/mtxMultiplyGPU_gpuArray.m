function [ C ] = mtxMultiplyGPU_gpuArray( A,B )
    [N,M] = size(A);
    [K,L] = size(B);
    assert(M == K,'Mtx multiply failed! Matrix dimensions do not match');
%     A = gpuArray(A);
%     B = gpuArray(B);
    C = A*B;

end


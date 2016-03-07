%% initialization
clear;clc;
N = 2048;
A = rand(N);
B = rand(N,N);

%% CPU version
tic
C_cpu = mtxMultiplyCPU(A,B);
time_cpu = toc;

%% GPU version
A = gpuArray(A);
B = gpuArray(B);
tic;
C_gpu = mtxMultiplyGPU_gpuArray(A,B);
time_gpuarray = toc;
speedup = time_cpu/time_gpuarray

%%
C_gpu = gather(C_gpu);
C = (C_cpu - C_gpu);
error = norm(C);
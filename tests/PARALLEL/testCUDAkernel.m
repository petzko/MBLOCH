K = 1024;
N = 1024;
L = 1024;

A = 2*ones(K,N);
B = 3*ones(N,L);

tic
C = A*B;
cpu_tm = toc;

A_g = 2*gpuArray(complex(ones(K,N)));
B_g = 3*gpuArray(complex(ones(N,L)));
C_g = gpuArray(complex(zeros(K,L)));

kernel = parallel.gpu.CUDAKernel('mtxmultiply.ptx','mtxmultiply.cu');
blockDim = [32 32];
gridDim  = [ceil(L/blockDim(1)) ceil(K/blockDim(2))];
kernel.GridSize = gridDim;
kernel.ThreadBlockSize = blockDim;

tic;
[C_g] = feval(kernel,A_g,B_g,C_g,N,K,L);
gpu_tm = toc;

speedup = cpu_tm/gpu_tm


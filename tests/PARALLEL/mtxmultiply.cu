__global__ void mtx_multiply(const double2 * A, const double2 * B, double2 * C,
	const int N, const int K, const int L ){


	int col = threadIdx.x + blockDim.x*blockIdx.x;
	int row = threadIdx.y + blockDim.y*blockIdx.y;

	if (col >= L || row >= K) 
		return;

	C[row*L + col].x = 0.0;
	C[row*L + col].y = 0.0;
	for(int k = 0; k < N ; k++){

		C[row*L + col].x += A[row*N+k].x*B[k*L+col].x - A[row*N+k].y*B[k*L+col].y;
		C[row*L + col].y += A[row*N+k].x*B[k*L+col].y + A[row*N+k].y*B[k*L+col].x;

	}
}
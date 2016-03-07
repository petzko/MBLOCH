#include <cuda.h>

__global__ void make_step(const double* D_N, const double2 * V, const double2* P,
							const double* c,
							const double dt, const int N, const int m, const int store_idx, 
							const int * order, 
							double2* F , double2* V_new)

{

	int t = threadIdx.x + blockDim.x*blockIdx.x;

	if(t >= N || t < 0)
		return;

	F[t*m + store_idx] = P[t];

	for(int k = 0; k < N; k++){
		F[t*m + store_idx].x += D_N[t*N+k]*V[k].x;
     	F[t*m + store_idx].y += D_N[t*N+k]*V[k].y;
	}

	V_new[t] = V[t];

	for(int s = 0; s < m; s++) {
		V_new[t].x += dt*c[order[s]]*F[t*m+s].x;	
		V_new[t].y += dt*c[order[s]]*F[t*m+s].y;	
	
	}

}
/* Implement Mixed Gradient (Lee et al. JCP 2016)*/
#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>
#include <stdio.h>

#define NBLOCKS 560
#define NTHREADS 128

void dvc_ScaLBL_D3Q19_MixedGradient(int *Map, double *Phi, double *Gradient, int start, int finish, int Np, int Nx, int Ny, int Nz,
                                    sycl::nd_item<3> item_ct1)
{
	// BUILD_FIX_AFTER_MIGRATION : Replaced static int with static const int
	static const int D3Q19[18][3]={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},
			{1,1,0},{-1,-1,0},{1,-1,0},{-1,1,0},
			{1,0,1},{-1,0,-1},{1,0,-1},{-1,0,1},
			{0,1,1},{0,-1,-1},{0,1,-1},{0,-1,1}};

	int i,j,k,n,N,idx;
	int np,np2,nm; // neighbors
	double v,vp,vp2,vm; // values at neighbors
	double grad;
	N = Nx*Ny*Nz;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){

		//........Get 1-D index for this thread....................
                idx = start +
                      S * item_ct1.get_group(2) * item_ct1.get_local_range(2) +
                      s * item_ct1.get_local_range(2) +
                      item_ct1.get_local_id(2);

                if (idx<finish){
			n = Map[idx]; // layout in regular array
			//.......Back out the 3-D indices for node n..............
			k = n/(Nx*Ny);
			j = (n-Nx*Ny*k)/Nx;
			i = n-Nx*Ny*k-Nx*j;	
			v = Phi[n];
			grad = 0.0;
			for (int q=0; q<6; q++){
				int iqx = D3Q19[q][0];
				int iqy = D3Q19[q][1];
				int iqz = D3Q19[q][2];
				np = (k+iqz)*Nx*Ny + (j+iqy)*Nx + i + iqx;
				np2 = (k+2*iqz)*Nx*Ny + (j+2*iqy)*Nx + i + 2*iqx;
				nm = (k-iqz)*Nx*Ny + (j-iqy)*Nx + i - iqx;
				vp = Phi[np];
				vp2 = Phi[np2];
				vm = Phi[nm];
				grad += 0.25*(5.0*vp-vp2-3.0*v-vm);
			}
			for (int q=6; q<18; q++){
				int iqx = D3Q19[q][0];
				int iqy = D3Q19[q][1];
				int iqz = D3Q19[q][2];
				np = (k+iqz)*Nx*Ny + (j+iqy)*Nx + i + iqx;
				np2 = (k+2*iqz)*Nx*Ny + (j+2*iqy)*Nx + i + 2*iqx;
				nm = (k-iqz)*Nx*Ny + (j-iqy)*Nx + i - iqx;
				vp = Phi[np];
				vp2 = Phi[np2];
				vm = Phi[nm];
				grad += 0.125*(5.0*vp-vp2-3.0*v-vm);
			}
			Gradient[n] = grad;
		}
	}
}

extern "C" void ScaLBL_D3Q19_MixedGradient(int *Map, double *Phi, double *Gradient, int start, int finish, int Np, int Nx, int Ny, int Nz)
{
        /*
        DPCT1007:82: Migration of cudaProfilerStart is not supported.
        */
		// BUILD_FIX_AFTER_MIGRATION : Commented below line.
        //cudaProfilerStart();
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, NBLOCKS) *
                                  sycl::range<3>(1, 1, NTHREADS),
                              sycl::range<3>(1, 1, NTHREADS)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q19_MixedGradient(Map, Phi, Gradient, start,
                                                   finish, Np, Nx, Ny, Nz,
                                                   item_ct1);
            });

        /*
        DPCT1010:84: SYCL uses exceptions to report errors and does not use the
        error codes. The call was replaced with 0. You need to rewrite this
        code.
        */
        int err = 0;

        /*
        DPCT1007:83: Migration of cudaProfilerStop is not supported.
        */
		// BUILD_FIX_AFTER_MIGRATION : Commented below line.
        //cudaProfilerStop();
}


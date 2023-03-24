// GPU Functions for D3Q7 Lattice Boltzmann Methods
#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>
#include <stdio.h>

#define NBLOCKS 560
#define NTHREADS 128

void dvc_ScaLBL_Scalar_Pack(int *list, int count, double *sendbuf, double *Data, int N,
                            sycl::nd_item<3> item_ct1){
	//....................................................................................
	// Pack distribution q into the send buffer for the listed lattice sites
	// dist may be even or odd distributions stored by stream layout
	//....................................................................................
	int idx,n;
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
        if (idx<count){
		n = list[idx];
		sendbuf[idx] = Data[n];
	}
}
void dvc_ScaLBL_Scalar_Unpack(int *list, int count, double *recvbuf, double *Data, int N,
                              sycl::nd_item<3> item_ct1){
	//....................................................................................
	// Pack distribution q into the send buffer for the listed lattice sites
	// dist may be even or odd distributions stored by stream layout
	//....................................................................................
	int idx,n;
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
        if (idx<count){
		n = list[idx];
		Data[n] = recvbuf[idx];
	}
}

void dvc_ScaLBL_PackDenD3Q7(int *list, int count, double *sendbuf, int number, double *Data, int N,
                            sycl::nd_item<3> item_ct1){
	//....................................................................................
	// Pack distribution into the send buffer for the listed lattice sites
	//....................................................................................
	int idx,n,component;
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
        if (idx<count){
		for (component=0; component<number; component++){
			n = list[idx];
			sendbuf[idx*number+component] = Data[number*n+component];
			Data[number*n+component] = 0.0;	// Set the data value to zero once it's in the buffer!
		}
	}
}


void dvc_ScaLBL_UnpackDenD3Q7(int *list, int count, double *recvbuf, int number, double *Data, int N,
                              sycl::nd_item<3> item_ct1){
	//....................................................................................
	// Unack distribution from the recv buffer
	// Sum to the existing density value
	//....................................................................................
	int idx,n,component;
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
        if (idx<count){
			for (component=0; component<number; component++){
			n = list[idx];
			Data[number*n+component] += recvbuf[idx*number+component];
		}
	}
}

void dvc_ScaLBL_D3Q7_Unpack(int q,  int *list,  int start, int count,
		double *recvbuf, double *dist, int N, sycl::nd_item<3> item_ct1){
	//....................................................................................
	// Unpack distribution from the recv buffer
	// Distribution q matche Cqx, Cqy, Cqz
	// swap rule means that the distributions in recvbuf are OPPOSITE of q
	// dist may be even or odd distributions stored by stream layout
	//....................................................................................
	int n,idx;
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
        if (idx<count){
		// Get the value from the list -- note that n is the index is from the send (non-local) process
		n = list[idx];
		// unpack the distribution to the proper location
		if (!(n<0)) { dist[q*N+n] = recvbuf[start+idx];
		//printf("%f \n",,dist[q*N+n]);
		}
	}
}

void dvc_ScaLBL_D3Q7_Reflection_BC_z(int *list, double *dist, int count, int Np,
                                     sycl::nd_item<3> item_ct1){
	int idx, n;
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
        if (idx < count){
		n = list[idx];
		double f5 = 0.222222222222222222222222 - dist[6*Np+n];
		dist[6*Np+n] = f5;
	}
}

void dvc_ScaLBL_D3Q7_Reflection_BC_Z(int *list, double *dist, int count, int Np,
                                     sycl::nd_item<3> item_ct1){
	int idx, n;
        idx = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
              item_ct1.get_local_id(2);
        if (idx < count){
		n = list[idx];
		double f6 = 0.222222222222222222222222 - dist[5*Np+n];
		dist[5*Np+n] = f6;
	}
}
void dvc_ScaLBL_D3Q7_Init(char *ID, double *f_even, double *f_odd, double *Den, int Nx, int Ny, int Nz,
                          sycl::nd_item<3> item_ct1)
{
	int n,N;
	N = Nx*Ny*Nz;
	double value;
	char id;
	int S = N/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
                n = S * item_ct1.get_group(2) * item_ct1.get_local_range(2) +
                    s * item_ct1.get_local_range(2) + item_ct1.get_local_id(2);
                if (n<N){
		   id = ID[n];
			if (id > 0){
				value = Den[n];
				f_even[n] = 0.3333333333333333*value;
				f_odd[n] = 0.1111111111111111*value;		//double(100*n)+1.f;
				f_even[N+n] = 0.1111111111111111*value;	//double(100*n)+2.f;
				f_odd[N+n] = 0.1111111111111111*value;	//double(100*n)+3.f;
				f_even[2*N+n] = 0.1111111111111111*value;	//double(100*n)+4.f;
				f_odd[2*N+n] = 0.1111111111111111*value;	//double(100*n)+5.f;
				f_even[3*N+n] = 0.1111111111111111*value;	//double(100*n)+6.f;
			}
			else{
				for(int q=0; q<3; q++){
					f_even[q*N+n] = -1.0;
					f_odd[q*N+n] = -1.0;
				}
				f_even[3*N+n] = -1.0;
			}
		}
	}
}

//*************************************************************************
void dvc_ScaLBL_D3Q7_Swap(char *ID, double *disteven, double *distodd, int Nx, int Ny, int Nz,
                          sycl::nd_item<3> item_ct1)
{
	int i,j,k,n,nn,N;
	// distributions
	double f1,f2,f3,f4,f5,f6;
	char id;
	N = Nx*Ny*Nz;

	int S = N/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
                n = S * item_ct1.get_group(2) * item_ct1.get_local_range(2) +
                    s * item_ct1.get_local_range(2) + item_ct1.get_local_id(2);

                if (n<N ){
		   	id = ID[n];
			if (id > 0){
                //.......Back out the 3-D indices for node n..............
                k = n/(Nx*Ny);
                j = (n-Nx*Ny*k)/Nx;
                i = n-Nx*Ny*k-Nx*j;
                //........................................................................
                // Retrieve even distributions from the local node (swap convention)
                //		f0 = disteven[n];  // Does not particupate in streaming
                f1 = distodd[n];
                f3 = distodd[N+n];
                f5 = distodd[2*N+n];
                //........................................................................

                //........................................................................
                // Retrieve odd distributions from neighboring nodes (swap convention)
                //........................................................................
                nn = n+1;							// neighbor index (pull convention)
                if (!(i+1<Nx))	nn -= Nx;			// periodic BC along the x-boundary
                //if (i+1<Nx){
                f2 = disteven[N+nn];					// pull neighbor for distribution 2
                if (!(f2 < 0.0)){
                    distodd[n] = f2;
                    disteven[N+nn] = f1;
                }
                //}
                //........................................................................
                nn = n+Nx;							// neighbor index (pull convention)
                if (!(j+1<Ny))	nn -= Nx*Ny;		// Perioidic BC along the y-boundary
                //if (j+1<Ny){
                f4 = disteven[2*N+nn];				// pull neighbor for distribution 4
                if (!(f4 < 0.0)){
                    distodd[N+n] = f4;
                    disteven[2*N+nn] = f3;
                }
                //........................................................................
                nn = n+Nx*Ny;						// neighbor index (pull convention)
                if (!(k+1<Nz))	nn -= Nx*Ny*Nz;		// Perioidic BC along the z-boundary
                //if (k+1<Nz){
                f6 = disteven[3*N+nn];				// pull neighbor for distribution 6
                if (!(f6 < 0.0)){
                    distodd[2*N+n] = f6;
                    disteven[3*N+nn] = f5;
                }
			}
		}
	}
}

//*************************************************************************
void dvc_ScaLBL_D3Q7_Density(char *ID, double *disteven, double *distodd, double *Den,
		int Nx, int Ny, int Nz, sycl::nd_item<3> item_ct1)
{
	char id;
	int n;
	double f0,f1,f2,f3,f4,f5,f6;
	int N = Nx*Ny*Nz;

	int S = N/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
                n = S * item_ct1.get_group(2) * item_ct1.get_local_range(2) +
                    s * item_ct1.get_local_range(2) + item_ct1.get_local_id(2);
                if (n<N){
            id = ID[n];
            if (id > 0 ){
                // Read the distributions
                f0 = disteven[n];
                f2 = disteven[N+n];
                f4 = disteven[2*N+n];
                f6 = disteven[3*N+n];
                f1 = distodd[n];
                f3 = distodd[N+n];
                f5 = distodd[2*N+n];
                // Compute the density
                Den[n] = f0+f1+f2+f3+f4+f5+f6;
            }
		}
	}
}

extern "C" void ScaLBL_D3Q7_Reflection_BC_z(int *list, double *dist, int count, int Np){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:75: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_Reflection_BC_z(list, dist, count, Np,
                                                    item_ct1);
            });
        /*
        DPCT1010:311: SYCL uses exceptions to report errors and does not use the
        error codes. The call was replaced with 0. You need to rewrite this
        code.
        */
        int err = 0;
}

extern "C" void ScaLBL_D3Q7_Reflection_BC_Z(int *list, double *dist, int count, int Np){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:76: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_Reflection_BC_Z(list, dist, count, Np,
                                                    item_ct1);
            });
        /*
        DPCT1010:313: SYCL uses exceptions to report errors and does not use the
        error codes. The call was replaced with 0. You need to rewrite this
        code.
        */
        int err = 0;
}

extern "C" void ScaLBL_D3Q7_Unpack(int q, int *list,  int start, int count, double *recvbuf, double *dist, int N){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:77: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_Unpack(q, list, start, count, recvbuf, dist,
                                           N, item_ct1);
            });
}

extern "C" void ScaLBL_Scalar_Pack(int *list, int count, double *sendbuf, double *Data, int N){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:78: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_Scalar_Pack(list, count, sendbuf, Data, N,
                                           item_ct1);
            });
}

extern "C" void ScaLBL_Scalar_Unpack(int *list, int count, double *recvbuf, double *Data, int N){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:79: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_Scalar_Unpack(list, count, recvbuf, Data, N,
                                             item_ct1);
            });
}
extern "C" void ScaLBL_PackDenD3Q7(int *list, int count, double *sendbuf, int number, double *Data, int N){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:80: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_PackDenD3Q7(list, count, sendbuf, number, Data,
                                           N, item_ct1);
            });
}

extern "C" void ScaLBL_UnpackDenD3Q7(int *list, int count, double *recvbuf, int number, double *Data, int N){
	int GRID = count / 512 + 1;
        /*
        DPCT1049:81: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, 512),
                              sycl::range<3>(1, 1, 512)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_UnpackDenD3Q7(list, count, recvbuf, number, Data,
                                             N, item_ct1);
            });
}

extern "C" void ScaLBL_D3Q7_Init(char *ID, double *f_even, double *f_odd, double *Den, int Nx, int Ny, int Nz){
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, NBLOCKS) *
                                  sycl::range<3>(1, 1, NTHREADS),
                              sycl::range<3>(1, 1, NTHREADS)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_Init(ID, f_even, f_odd, Den, Nx, Ny, Nz,
                                         item_ct1);
            });
}

extern "C" void ScaLBL_D3Q7_Swap(char *ID, double *disteven, double *distodd, int Nx, int Ny, int Nz){
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, NBLOCKS) *
                                  sycl::range<3>(1, 1, NTHREADS),
                              sycl::range<3>(1, 1, NTHREADS)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_Swap(ID, disteven, distodd, Nx, Ny, Nz,
                                         item_ct1);
            });
}

extern "C" void ScaLBL_D3Q7_Density(char *ID, double *disteven, double *distodd, double *Den,
										int Nx, int Ny, int Nz){
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, NBLOCKS) *
                                  sycl::range<3>(1, 1, NTHREADS),
                              sycl::range<3>(1, 1, NTHREADS)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_Density(ID, disteven, distodd, Den, Nx, Ny,
                                            Nz, item_ct1);
            });
}


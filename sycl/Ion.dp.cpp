#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>
#include <stdio.h>
#include <math.h>
//#include <cuda_profiler_api.h>

#define NBLOCKS 1024
#define NTHREADS 512

extern "C" void Membrane_D3Q19_Unpack(int q, int *list, int *links, int start, int linkCount,
                                    double *recvbuf, double *dist, int N) {
    //....................................................................................
    // Unack distribution from the recv buffer
    // Distribution q matche Cqx, Cqy, Cqz
    // swap rule means that the distributions in recvbuf are OPPOSITE of q
    // dist may be even or odd distributions stored by stream layout
    //....................................................................................
    int n, idx, link;
    for (link=0; link<linkCount; link++){

    	idx = links[start+link];
        // Get the value from the list -- note that n is the index is from the send (non-local) process
        n = list[start + idx];
        // unpack the distribution to the proper location
        if (!(n < 0))
            dist[q * N + n] = recvbuf[start + idx];
    }
}

extern "C" void Membrane_D3Q19_Transport(int q, int *list, int *links, double *coef, int start, int offset, 
		int linkCount, double *recvbuf, double *dist, int N){
    //....................................................................................
    // Unack distribution from the recv buffer
    // Distribution q matche Cqx, Cqy, Cqz
    // swap rule means that the distributions in recvbuf are OPPOSITE of q
    // dist may be even or odd distributions stored by stream layout
    //....................................................................................
    int n, idx, link;
    double alpha;
    for (link=offset; link<linkCount; link++){

    	idx = list[start+link];
        // Get the value from the list -- note that n is the index is from the send (non-local) process
        n = list[start + idx];
        alpha = coef[start + idx];
        // unpack the distribution to the proper location
        if (!(n < 0))
            dist[q * N + n] = alpha*recvbuf[start + idx];
    }
}

void dvc_ScaLBL_D3Q7_Membrane_AssignLinkCoef(int *membrane, int *Map, double *Distance, double *Psi, double *coef,
		double Threshold, double MassFractionIn, double MassFractionOut, double ThresholdMassFractionIn, double ThresholdMassFractionOut,
		int memLinks, int Nx, int Ny, int Nz, int Np, sycl::nd_item<3> item_ct1){

	int link,iq,ip,nq,np,nqm,npm;
	double aq, ap, membranePotential;
	/* Interior Links */

	int S = memLinks/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
                link = S * item_ct1.get_group(2) * item_ct1.get_local_range(2) +
                       s * item_ct1.get_local_range(2) +
                       item_ct1.get_local_id(2);
                if (link < memLinks) {

			// inside             	//outside
			aq = MassFractionIn;	ap = MassFractionOut;  
			iq = membrane[2*link]; 	ip = membrane[2*link+1];
			nq = iq%Np;				np = ip%Np;
			nqm = Map[nq];			npm = Map[np]; // strided layout

			/* membrane potential for this link */
			membranePotential = Psi[nqm] - Psi[npm];
			if (membranePotential > Threshold){
				aq = ThresholdMassFractionIn;	ap = ThresholdMassFractionOut;  
			}

			/* Save the mass transfer coefficients */
			coef[2*link] = aq;		coef[2*link+1] = ap;
		}
	}
}

void dvc_ScaLBL_D3Q7_Membrane_AssignLinkCoef_halo(
		const int Cqx, const int Cqy, int const Cqz, 
		int *Map, double *Distance, double *Psi, double Threshold, 
		double MassFractionIn, double MassFractionOut, double ThresholdMassFractionIn, double ThresholdMassFractionOut,
		int *d3q7_recvlist, int *d3q7_linkList, double *coef, int start, int nlinks, int count,
		const int N, const int Nx, const int Ny, const int Nz,
		sycl::nd_item<3> item_ct1) {
	//....................................................................................
	// Unack distribution from the recv buffer
	// Distribution q matche Cqx, Cqy, Cqz
	// swap rule means that the distributions in recvbuf are OPPOSITE of q
	// dist may be even or odd distributions stored by stream layout
	//....................................................................................
	int n, idx, nqm, npm, label, i, j, k;
	double distanceLocal, distanceNonlocal;
	double psiLocal, psiNonlocal, membranePotential;
	double ap,aq; // coefficient

	/* second enforce custom rule for membrane links */
	int S = (count-nlinks)/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
                idx = S * item_ct1.get_group(2) * item_ct1.get_local_range(2) +
                      s * item_ct1.get_local_range(2) +
                      item_ct1.get_local_id(2);

                if (idx < count) {

			n = d3q7_recvlist[idx];
			label = d3q7_linkList[idx];
			ap = 1.0;  // regular streaming rule
			aq = 1.0;
			if (label > 0 && !(n < 0)){
				nqm = Map[n];
				distanceLocal = Distance[nqm];  
				psiLocal = Psi[nqm];

				// Get the 3-D indices from the send process
				k = nqm/(Nx*Ny); j = (nqm-Nx*Ny*k)/Nx; i = nqm-Nx*Ny*k-Nx*j;
				// Streaming link the non-local distribution
				i -= Cqx; j -= Cqy; k -= Cqz;
				npm = k*Nx*Ny + j*Nx + i;
				distanceNonlocal = Distance[npm];  
				psiNonlocal = Psi[npm];

				membranePotential = psiLocal - psiNonlocal;
				aq = MassFractionIn;
				ap = MassFractionOut;

				/* link is inside membrane */
				if (distanceLocal > 0.0){
					if (membranePotential < Threshold*(-1.0)){
						ap = MassFractionIn;
						aq = MassFractionOut;
					}
					else {
						ap = ThresholdMassFractionIn;
						aq = ThresholdMassFractionOut;
					}
				}
				else if (membranePotential > Threshold){
					aq = ThresholdMassFractionIn;
					ap = ThresholdMassFractionOut;
				}
			}
			coef[2*idx]=aq;
			coef[2*idx+1]=ap;
		}
	}
}

void dvc_ScaLBL_D3Q7_Membrane_Unpack(int q,  
		int *d3q7_recvlist, double *recvbuf, int count,
		double *dist, int N,  double *coef, sycl::nd_item<3> item_ct1)  {
	//....................................................................................
	// Unack distribution from the recv buffer
	// Distribution q matche Cqx, Cqy, Cqz
	// swap rule means that the distributions in recvbuf are OPPOSITE of q
	// dist may be even or odd distributions stored by stream layout
	//....................................................................................
	int n, idx, link;
	double fq,fp,fqq,ap,aq; // coefficient

	/* second enforce custom rule for membrane links */
	int S = count/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
                idx = S * item_ct1.get_group(2) * item_ct1.get_local_range(2) +
                      s * item_ct1.get_local_range(2) +
                      item_ct1.get_local_id(2);
                if (idx < count){
	    	n = d3q7_recvlist[idx];
	        // update link based on mass transfer coefficients
	        if (!(n < 0)){
	        	aq = coef[2*idx];
	        	ap = coef[2*idx+1];
	        	fq = dist[q * N + n];
	        	fp = recvbuf[idx];
	        	fqq = (1-aq)*fq+ap*fp;
	            dist[q * N + n] = fqq;
	        }
		} 
	}
}

void dvc_ScaLBL_D3Q7_Membrane_IonTransport(int *membrane, double *coef, 
		double *dist, double *Den, int memLinks, int Np, sycl::nd_item<3> item_ct1){	
	int link,iq,ip,nq,np;
	double aq, ap, fq, fp, fqq, fpp, Cq, Cp;

	int S = memLinks/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
                link = S * item_ct1.get_group(2) * item_ct1.get_local_range(2) +
                       s * item_ct1.get_local_range(2) +
                       item_ct1.get_local_id(2);
                if (link < memLinks){

			// inside             	//outside
			aq = coef[2*link];		ap = coef[2*link+1];
			iq = membrane[2*link]; 	ip = membrane[2*link+1];
			nq = iq%Np;				np = ip%Np;
			fq  = dist[iq];			fp = dist[ip];
			fqq = (1-aq)*fq+ap*fp;	fpp = (1-ap)*fp+aq*fq;
			Cq = Den[nq];			Cp = Den[np];
			Cq += fqq - fq;			Cp += fpp - fp;
			Den[nq] = Cq;			Den[np] = Cp;
			dist[iq] = fqq;			dist[ip] = fpp;
		}
	}
}


void dvc_ScaLBL_D3Q7_AAodd_IonConcentration(int *neighborList, double *dist, double *Den, int start, int finish, int Np,
                                            sycl::nd_item<3> item_ct1){
    int n,nread;
    double fq,Ci;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
                n = S * item_ct1.get_group(2) * item_ct1.get_local_range(2) +
                    s * item_ct1.get_local_range(2) + item_ct1.get_local_id(2) +
                    start;
                if (n<finish) {

            // q=0
            fq = dist[n];
            Ci = fq;

            // q=1
            nread = neighborList[n]; 
            fq = dist[nread]; 
            Ci += fq;
            
            // q=2
            nread = neighborList[n+Np]; 
            fq = dist[nread];  
            Ci += fq;
            
            // q=3
            nread = neighborList[n+2*Np]; 
            fq = dist[nread];
            Ci += fq;
            
            // q=4
            nread = neighborList[n+3*Np]; 
            fq = dist[nread];
            Ci += fq;
            
            // q=5
            nread = neighborList[n+4*Np];
            fq = dist[nread];
            Ci += fq;
            
            // q=6
            nread = neighborList[n+5*Np];
            fq = dist[nread];
            Ci += fq;

            Den[n]=Ci;
		}
	}
}

void dvc_ScaLBL_D3Q7_AAeven_IonConcentration(double *dist, double *Den, int start, int finish, int Np,
                                             sycl::nd_item<3> item_ct1){
    int n;
    double fq,Ci;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
                n = S * item_ct1.get_group(2) * item_ct1.get_local_range(2) +
                    s * item_ct1.get_local_range(2) + item_ct1.get_local_id(2) +
                    start;
                if (n<finish) {

            // q=0
            fq = dist[n];
            Ci = fq;
            
            // q=1
            fq = dist[2*Np+n];
            Ci += fq;

            // q=2
            fq = dist[1*Np+n];
            Ci += fq;

            // q=3
            fq = dist[4*Np+n];
            Ci += fq;

            // q=4
            fq = dist[3*Np+n];
            Ci += fq;

            // q=5
            fq = dist[6*Np+n];
            Ci += fq;

            // q=6
            fq = dist[5*Np+n];
            Ci += fq;

            Den[n]=Ci;
		}
	}
}

void dvc_ScaLBL_D3Q7_AAodd_Ion(int *neighborList, double *dist, double *Den, double *FluxDiffusive, double *FluxAdvective, double *FluxElectrical, double *Velocity, double *ElectricField, 
                                           double Di, int zi, double rlx, double Vt, int start, int finish, int Np,
                                           sycl::nd_item<3> item_ct1){
	int n;
	double Ci;
    double ux,uy,uz;
    double uEPx,uEPy,uEPz;//electrochemical induced velocity
    double Ex,Ey,Ez;//electrical field
    double flux_diffusive_x,flux_diffusive_y,flux_diffusive_z;
	double f0,f1,f2,f3,f4,f5,f6;
	double X,Y,Z,factor_x,factor_y,factor_z;
	int nr1,nr2,nr3,nr4,nr5,nr6;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
                n = S * item_ct1.get_group(2) * item_ct1.get_local_range(2) +
                    s * item_ct1.get_local_range(2) + item_ct1.get_local_id(2) +
                    start;
                if (n<finish) {

	        //Load data
	        Ex = ElectricField[n + 0 * Np];
	        Ey = ElectricField[n + 1 * Np];
	        Ez = ElectricField[n + 2 * Np];
	        ux = Velocity[n + 0 * Np];
	        uy = Velocity[n + 1 * Np];
	        uz = Velocity[n + 2 * Np];
	        uEPx = zi * Di / Vt * Ex;
	        uEPy = zi * Di / Vt * Ey;
	        uEPz = zi * Di / Vt * Ez;

	        // q=0
	        f0 = dist[n];
	        // q=1
	        nr1 = neighborList[n]; // neighbor 2 ( > 10Np => odd part of dist)
	        f1 = dist[nr1];        // reading the f1 data into register fq
	        // q=2
	        nr2 = neighborList[n + Np]; // neighbor 1 ( < 10Np => even part of dist)
	        f2 = dist[nr2];             // reading the f2 data into register fq
	        // q=3
	        nr3 = neighborList[n + 2 * Np]; // neighbor 4
	        f3 = dist[nr3];
	        // q=4
	        nr4 = neighborList[n + 3 * Np]; // neighbor 3
	        f4 = dist[nr4];
	        // q=5
	        nr5 = neighborList[n + 4 * Np];
	        f5 = dist[nr5];
	        // q=6
	        nr6 = neighborList[n + 5 * Np];
	        f6 = dist[nr6];

	        // compute diffusive flux
	        Ci = f0 + f1 + f2 + f3 + f4 + f5 + f6;
	        flux_diffusive_x = (1.0 - 0.5 * rlx) * ((f1 - f2) - ux * Ci);
	        flux_diffusive_y = (1.0 - 0.5 * rlx) * ((f3 - f4) - uy * Ci);
	        flux_diffusive_z = (1.0 - 0.5 * rlx) * ((f5 - f6) - uz * Ci);
	        FluxDiffusive[n + 0 * Np] = flux_diffusive_x;
	        FluxDiffusive[n + 1 * Np] = flux_diffusive_y;
	        FluxDiffusive[n + 2 * Np] = flux_diffusive_z;
	        FluxAdvective[n + 0 * Np] = ux * Ci;
	        FluxAdvective[n + 1 * Np] = uy * Ci;
	        FluxAdvective[n + 2 * Np] = uz * Ci;
	        FluxElectrical[n + 0 * Np] = uEPx * Ci;
	        FluxElectrical[n + 1 * Np] = uEPy * Ci;
	        FluxElectrical[n + 2 * Np] = uEPz * Ci;
	        
	        Den[n] = Ci;

	        /* use logistic function to prevent negative distributions*/
	        X = 4.0 * (ux + uEPx);
	        Y = 4.0 * (uy + uEPy);
	        Z = 4.0 * (uz + uEPz);
                factor_x = X / sycl::sqrt(1 + X * X);
                factor_y = Y / sycl::sqrt(1 + Y * Y);
                factor_z = Z / sycl::sqrt(1 + Z * Z);

                // q=0
	        dist[n] = f0 * (1.0 - rlx) + rlx * 0.25 * Ci;

	        // q = 1
	        dist[nr2] =
	            f1 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 + factor_x);
	        //f1 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 + 4.0 * (ux + uEPx));


	        // q=2
	        dist[nr1] =
	                f2 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 - factor_x);
	        //f2 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 - 4.0 * (ux + uEPx));

	        // q = 3
	        dist[nr4] =
	                f3 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 + factor_y );
	        //f3 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 + 4.0 * (uy + uEPy));

	        // q = 4
	        dist[nr3] =
	                f4 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 - factor_y);
	        //f4 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 - 4.0 * (uy + uEPy));

	        // q = 5
	        dist[nr6] =
	                f5 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 +  factor_z);
	        //f5 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 + 4.0 * (uz + uEPz));

	        // q = 6
	        dist[nr5] =
	            f6 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 - factor_z);

		}
	}
}

void dvc_ScaLBL_D3Q7_AAeven_Ion(double *dist, double *Den, double *FluxDiffusive, double *FluxAdvective, double *FluxElectrical, double *Velocity, double *ElectricField, 
                                            double Di, int zi, double rlx, double Vt, int start, int finish, int Np,
                                            sycl::nd_item<3> item_ct1){
	int n;
	double Ci;
    double ux,uy,uz;
    double uEPx,uEPy,uEPz;//electrochemical induced velocity
    double Ex,Ey,Ez;//electrical field
    double flux_diffusive_x,flux_diffusive_y,flux_diffusive_z;
	double f0,f1,f2,f3,f4,f5,f6;
	double X,Y,Z,factor_x,factor_y,factor_z;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
                n = S * item_ct1.get_group(2) * item_ct1.get_local_range(2) +
                    s * item_ct1.get_local_range(2) + item_ct1.get_local_id(2) +
                    start;
                if (n<finish) {

		       //Load data
		        //Ci = Den[n];
		        Ex = ElectricField[n + 0 * Np];
		        Ey = ElectricField[n + 1 * Np];
		        Ez = ElectricField[n + 2 * Np];
		        ux = Velocity[n + 0 * Np];
		        uy = Velocity[n + 1 * Np];
		        uz = Velocity[n + 2 * Np];
		        uEPx = zi * Di / Vt * Ex;
		        uEPy = zi * Di / Vt * Ey;
		        uEPz = zi * Di / Vt * Ez;

		        f0 = dist[n];
		        f1 = dist[2 * Np + n];
		        f2 = dist[1 * Np + n];
		        f3 = dist[4 * Np + n];
		        f4 = dist[3 * Np + n];
		        f5 = dist[6 * Np + n];
		        f6 = dist[5 * Np + n];

		        // compute diffusive flux
		        Ci = f0 + f1 + f2 + f3 + f4 + f5 + f6;
		        flux_diffusive_x = (1.0 - 0.5 * rlx) * ((f1 - f2) - ux * Ci);
		        flux_diffusive_y = (1.0 - 0.5 * rlx) * ((f3 - f4) - uy * Ci);
		        flux_diffusive_z = (1.0 - 0.5 * rlx) * ((f5 - f6) - uz * Ci);
		        FluxDiffusive[n + 0 * Np] = flux_diffusive_x;
		        FluxDiffusive[n + 1 * Np] = flux_diffusive_y;
		        FluxDiffusive[n + 2 * Np] = flux_diffusive_z;
		        FluxAdvective[n + 0 * Np] = ux * Ci;
		        FluxAdvective[n + 1 * Np] = uy * Ci;
		        FluxAdvective[n + 2 * Np] = uz * Ci;
		        FluxElectrical[n + 0 * Np] = uEPx * Ci;
		        FluxElectrical[n + 1 * Np] = uEPy * Ci;
		        FluxElectrical[n + 2 * Np] = uEPz * Ci;
		        
		        Den[n] = Ci;
		        
		        /* use logistic function to prevent negative distributions*/
		        X = 4.0 * (ux + uEPx);
		        Y = 4.0 * (uy + uEPy);
		        Z = 4.0 * (uz + uEPz);
                        factor_x = X / sycl::sqrt(1 + X * X);
                        factor_y = Y / sycl::sqrt(1 + Y * Y);
                        factor_z = Z / sycl::sqrt(1 + Z * Z);

                        // q=0
		        dist[n] = f0 * (1.0 - rlx) + rlx * 0.25 * Ci;

		        // q = 1
		        dist[1 * Np + n] =
		                f1 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 + factor_x);
		        //f1 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 + 4.0 * (ux + uEPx));

		        // q=2
		        dist[2 * Np + n] =
		                f2 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 - factor_x);
		        //f2 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 - 4.0 * (ux + uEPx));

		        // q = 3
		        dist[3 * Np + n] =
		                f3 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 + factor_y);
		        //f3 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 + 4.0 * (uy + uEPy));

		        // q = 4
		        dist[4 * Np + n] =
		                f4 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 - factor_y);
		        //f4 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 - 4.0 * (uy + uEPy));

		        // q = 5
		        dist[5 * Np + n] =
		                f5 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 + factor_z);
		        //f5 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 + 4.0 * (uz + uEPz));

		        // q = 6
		        dist[6 * Np + n] =
		                f6 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 - factor_z);
		        //f6 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 - 4.0 * (uz + uEPz));
		}
	}
}

void dvc_ScaLBL_D3Q7_Ion_Init(double *dist, double *Den, double DenInit, int Np,
                              sycl::nd_item<3> item_ct1){

	int n;
	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
                n = S * item_ct1.get_group(2) * item_ct1.get_local_range(2) +
                    s * item_ct1.get_local_range(2) + item_ct1.get_local_id(2);
                if (n<Np) {
            dist[0*Np+n] = 0.25*DenInit;
            dist[1*Np+n] = 0.125*DenInit;		
            dist[2*Np+n] = 0.125*DenInit;	
            dist[3*Np+n] = 0.125*DenInit;	
            dist[4*Np+n] = 0.125*DenInit;	
            dist[5*Np+n] = 0.125*DenInit;	
            dist[6*Np+n] = 0.125*DenInit;	
            Den[n] = DenInit;
		}
	}
}

void dvc_ScaLBL_D3Q7_Ion_Init_FromFile(double *dist, double *Den, int Np,
                                       sycl::nd_item<3> item_ct1){

	int n;
    double DenInit;
	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
                n = S * item_ct1.get_group(2) * item_ct1.get_local_range(2) +
                    s * item_ct1.get_local_range(2) + item_ct1.get_local_id(2);
                if (n<Np) {
            DenInit = Den[n];
            dist[0*Np+n] = 0.25*DenInit;
            dist[1*Np+n] = 0.125*DenInit;		
            dist[2*Np+n] = 0.125*DenInit;	
            dist[3*Np+n] = 0.125*DenInit;	
            dist[4*Np+n] = 0.125*DenInit;	
            dist[5*Np+n] = 0.125*DenInit;	
            dist[6*Np+n] = 0.125*DenInit;	
		}
	}
}

void dvc_ScaLBL_D3Q7_Ion_ChargeDensity(double *Den, double *ChargeDensity, double IonValence, int ion_component, int start, int finish, int Np,
                                       sycl::nd_item<3> item_ct1){

    int n;
    double Ci;//ion concentration of species i
    double CD;//charge density
    double CD_tmp;
    double F = 96485.0;//Faraday's constant; unit[C/mol]; F=e*Na, where Na is the Avogadro constant

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
                n = S * item_ct1.get_group(2) * item_ct1.get_local_range(2) +
                    s * item_ct1.get_local_range(2) + item_ct1.get_local_id(2) +
                    start;
                if (n<finish) {

            Ci = Den[n+ion_component*Np];
            CD = ChargeDensity[n];
            if (ion_component == 0) CD=0.0;
            CD_tmp = F*IonValence*Ci;
            ChargeDensity[n] = CD + CD_tmp;

 //           Ci = Den[n+ion_component*Np];
   //         CD = ChargeDensity[n];
     //       CD_tmp = F*IonValence*Ci;
       //     ChargeDensity[n] = CD*(ion_component>0) + CD_tmp;
		}
	}
}
void dvc_ScaLBL_D3Q7_AAodd_Ion_v0(int *neighborList, double *dist,
                                      double *Den, double *FluxDiffusive,
                                      double *FluxAdvective,
                                      double *FluxElectrical, double *Velocity,
                                      double *ElectricField, double Di, int zi,
                                      double rlx, double Vt, int start,
                                      int finish, int Np,
                                      sycl::nd_item<3> item_ct1) {
    int n;
    double Ci;
    double ux, uy, uz;
    double uEPx, uEPy, uEPz; //electrochemical induced velocity
    double Ex, Ey, Ez;       //electrical field
    double flux_diffusive_x, flux_diffusive_y, flux_diffusive_z;
    double f0, f1, f2, f3, f4, f5, f6;
    //double X,Y,Z,factor_x, factor_y, factor_z;
    int nr1, nr2, nr3, nr4, nr5, nr6;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
                n = S * item_ct1.get_group(2) * item_ct1.get_local_range(2) +
                    s * item_ct1.get_local_range(2) + item_ct1.get_local_id(2) +
                    start;
                if (n<finish) {


        //Load data
        Ci = Den[n];
        Ex = ElectricField[n + 0 * Np];
        Ey = ElectricField[n + 1 * Np];
        Ez = ElectricField[n + 2 * Np];
        ux = Velocity[n + 0 * Np];
        uy = Velocity[n + 1 * Np];
        uz = Velocity[n + 2 * Np];
        uEPx = zi * Di / Vt * Ex;
        uEPy = zi * Di / Vt * Ey;
        uEPz = zi * Di / Vt * Ez;

        // q=0
        f0 = dist[n];
        // q=1
        nr1 = neighborList[n]; // neighbor 2 ( > 10Np => odd part of dist)
        f1 = dist[nr1];        // reading the f1 data into register fq
        // q=2
        nr2 = neighborList[n + Np]; // neighbor 1 ( < 10Np => even part of dist)
        f2 = dist[nr2];             // reading the f2 data into register fq
        // q=3
        nr3 = neighborList[n + 2 * Np]; // neighbor 4
        f3 = dist[nr3];
        // q=4
        nr4 = neighborList[n + 3 * Np]; // neighbor 3
        f4 = dist[nr4];
        // q=5
        nr5 = neighborList[n + 4 * Np];
        f5 = dist[nr5];
        // q=6
        nr6 = neighborList[n + 5 * Np];
        f6 = dist[nr6];

        // compute diffusive flux
        //Ci = f0 + f1 + f2 + f3 + f4 + f5 + f6;
        flux_diffusive_x = (1.0 - 0.5 * rlx) * ((f1 - f2) - ux * Ci);
        flux_diffusive_y = (1.0 - 0.5 * rlx) * ((f3 - f4) - uy * Ci);
        flux_diffusive_z = (1.0 - 0.5 * rlx) * ((f5 - f6) - uz * Ci);
        FluxDiffusive[n + 0 * Np] = flux_diffusive_x;
        FluxDiffusive[n + 1 * Np] = flux_diffusive_y;
        FluxDiffusive[n + 2 * Np] = flux_diffusive_z;
        FluxAdvective[n + 0 * Np] = ux * Ci;
        FluxAdvective[n + 1 * Np] = uy * Ci;
        FluxAdvective[n + 2 * Np] = uz * Ci;
        FluxElectrical[n + 0 * Np] = uEPx * Ci;
        FluxElectrical[n + 1 * Np] = uEPy * Ci;
        FluxElectrical[n + 2 * Np] = uEPz * Ci;
        
        //Den[n] = Ci;

        /* use logistic function to prevent negative distributions*/
        //X = 4.0 * (ux + uEPx);
        //Y = 4.0 * (uy + uEPy);
        //Z = 4.0 * (uz + uEPz);
        //factor_x = X / sqrt(1 + X*X);
        //factor_y = Y / sqrt(1 + Y*Y);
        //factor_z = Z / sqrt(1 + Z*Z);

        // q=0
        dist[n] = f0 * (1.0 - rlx) + rlx * 0.25 * Ci;

        // q = 1
        dist[nr2] =
        f1 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 + 4.0 * (ux + uEPx));
        //    f1 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 + factor_x);


        // q=2
        dist[nr1] =
        f2 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 - 4.0 * (ux + uEPx));
        //        f2 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 - factor_x);

        // q = 3
        dist[nr4] =
        f3 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 + 4.0 * (uy + uEPy));
        //        f3 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 + factor_y );

        // q = 4
        dist[nr3] =
        f4 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 - 4.0 * (uy + uEPy));
        //        f4 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 - factor_y);

        // q = 5
        dist[nr6] =
        f5 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 + 4.0 * (uz + uEPz));
        //        f5 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 +  factor_z);

        // q = 6
        dist[nr5] =
        f6 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 - 4.0 * (uz + uEPz));
        //    f6 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 - factor_z);

		}
    }
}

void dvc_ScaLBL_D3Q7_AAeven_Ion_v0(
    double *dist, double *Den, double *FluxDiffusive, double *FluxAdvective,
    double *FluxElectrical, double *Velocity, double *ElectricField, double Di,
    int zi, double rlx, double Vt, int start, int finish, int Np,
    sycl::nd_item<3> item_ct1) {
    int n;
    double Ci;
    double ux, uy, uz;
    double uEPx, uEPy, uEPz; //electrochemical induced velocity
    double Ex, Ey, Ez;       //electrical field
    double flux_diffusive_x, flux_diffusive_y, flux_diffusive_z;
    double f0, f1, f2, f3, f4, f5, f6;
    //double X,Y,Z, factor_x, factor_y, factor_z;

	int S = Np/NBLOCKS/NTHREADS + 1;
	for (int s=0; s<S; s++){
		//........Get 1-D index for this thread....................
                n = S * item_ct1.get_group(2) * item_ct1.get_local_range(2) +
                    s * item_ct1.get_local_range(2) + item_ct1.get_local_id(2) +
                    start;
                if (n<finish) {

        //Load data
        Ci = Den[n];
        Ex = ElectricField[n + 0 * Np];
        Ey = ElectricField[n + 1 * Np];
        Ez = ElectricField[n + 2 * Np];
        ux = Velocity[n + 0 * Np];
        uy = Velocity[n + 1 * Np];
        uz = Velocity[n + 2 * Np];
        uEPx = zi * Di / Vt * Ex;
        uEPy = zi * Di / Vt * Ey;
        uEPz = zi * Di / Vt * Ez;

        f0 = dist[n];
        f1 = dist[2 * Np + n];
        f2 = dist[1 * Np + n];
        f3 = dist[4 * Np + n];
        f4 = dist[3 * Np + n];
        f5 = dist[6 * Np + n];
        f6 = dist[5 * Np + n];

        // compute diffusive flux
        //Ci = f0 + f1 + f2 + f3 + f4 + f5 + f6;
        flux_diffusive_x = (1.0 - 0.5 * rlx) * ((f1 - f2) - ux * Ci);
        flux_diffusive_y = (1.0 - 0.5 * rlx) * ((f3 - f4) - uy * Ci);
        flux_diffusive_z = (1.0 - 0.5 * rlx) * ((f5 - f6) - uz * Ci);
        FluxDiffusive[n + 0 * Np] = flux_diffusive_x;
        FluxDiffusive[n + 1 * Np] = flux_diffusive_y;
        FluxDiffusive[n + 2 * Np] = flux_diffusive_z;
        FluxAdvective[n + 0 * Np] = ux * Ci;
        FluxAdvective[n + 1 * Np] = uy * Ci;
        FluxAdvective[n + 2 * Np] = uz * Ci;
        FluxElectrical[n + 0 * Np] = uEPx * Ci;
        FluxElectrical[n + 1 * Np] = uEPy * Ci;
        FluxElectrical[n + 2 * Np] = uEPz * Ci;
        
        //Den[n] = Ci;
        
        /* use logistic function to prevent negative distributions*/
        //X = 4.0 * (ux + uEPx);
        //Y = 4.0 * (uy + uEPy);
        //Z = 4.0 * (uz + uEPz);
        //factor_x = X / sqrt(1 + X*X);
        //factor_y = Y / sqrt(1 + Y*Y);
        //factor_z = Z / sqrt(1 + Z*Z);

        // q=0
        dist[n] = f0 * (1.0 - rlx) + rlx * 0.25 * Ci;

        // q = 1
        dist[1 * Np + n] =
        f1 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 + 4.0 * (ux + uEPx));
        //        f1 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 + factor_x);

        // q=2
        dist[2 * Np + n] =
        f2 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 - 4.0 * (ux + uEPx));
        //        f2 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 - factor_x);

        // q = 3
        dist[3 * Np + n] =
        f3 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 + 4.0 * (uy + uEPy));
        //        f3 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 + factor_y);

        // q = 4
        dist[4 * Np + n] =
        f4 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 - 4.0 * (uy + uEPy));
        //        f4 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 - factor_y);

        // q = 5
        dist[5 * Np + n] =
        f5 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 + 4.0 * (uz + uEPz));
        //        f5 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 + factor_z);

        // q = 6
        dist[6 * Np + n] =
        f6 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 - 4.0 * (uz + uEPz));
        //        f6 * (1.0 - rlx) + rlx * 0.125 * Ci * (1.0 - factor_z);
        }
    }
}

extern "C" void ScaLBL_D3Q7_AAeven_Ion_v0(
    double *dist, double *Den, double *FluxDiffusive, double *FluxAdvective,
    double *FluxElectrical, double *Velocity, double *ElectricField, double Di,
    int zi, double rlx, double Vt, int start, int finish, int Np) {

    /*
    DPCT1049:0: The work-group size passed to the SYCL kernel may exceed the
    limit. To get the device limit, query info::device::max_work_group_size.
    Adjust the work-group size if needed.
    */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, NBLOCKS) *
                                  sycl::range<3>(1, 1, NTHREADS),
                              sycl::range<3>(1, 1, NTHREADS)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_AAeven_Ion_v0(
                        dist, Den, FluxDiffusive, FluxAdvective, FluxElectrical,
                        Velocity, ElectricField, Di, zi, rlx, Vt, start, finish,
                        Np, item_ct1);
            });

        /*
        DPCT1010:100: SYCL uses exceptions to report errors and does not use the
        error codes. The call was replaced with 0. You need to rewrite this
        code.
        */
        int err = 0;
}

extern "C" void ScaLBL_D3Q7_AAodd_Ion_v0(int *neighborList, double *dist,
                                      double *Den, double *FluxDiffusive,
                                      double *FluxAdvective,
                                      double *FluxElectrical, double *Velocity,
                                      double *ElectricField, double Di, int zi,
                                      double rlx, double Vt, int start,
                                      int finish, int Np) {

        /*
        DPCT1049:1: The work-group size passed to the SYCL kernel may exceed the
        limit. To get the device limit, query info::device::max_work_group_size.
        Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, NBLOCKS) *
                                  sycl::range<3>(1, 1, NTHREADS),
                              sycl::range<3>(1, 1, NTHREADS)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_AAodd_Ion_v0(
                        neighborList, dist, Den, FluxDiffusive, FluxAdvective,
                        FluxElectrical, Velocity, ElectricField, Di, zi, rlx,
                        Vt, start, finish, Np, item_ct1);
            });

        /*
        DPCT1010:102: SYCL uses exceptions to report errors and does not use the
        error codes. The call was replaced with 0. You need to rewrite this
        code.
        */
        int err = 0;
} 
                                 


extern "C" void ScaLBL_D3Q7_AAodd_IonConcentration(int *neighborList, double *dist, double *Den, int start, int finish, int Np){

	//cudaProfilerStart();
        /*
        DPCT1049:2: The work-group size passed to the SYCL kernel may exceed the
        limit. To get the device limit, query info::device::max_work_group_size.
        Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, NBLOCKS) *
                                  sycl::range<3>(1, 1, NTHREADS),
                              sycl::range<3>(1, 1, NTHREADS)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_AAodd_IonConcentration(
                        neighborList, dist, Den, start, finish, Np, item_ct1);
            });

        /*
        DPCT1010:104: SYCL uses exceptions to report errors and does not use the
        error codes. The call was replaced with 0. You need to rewrite this
        code.
        */
        int err = 0;

        //cudaProfilerStop();
}

extern "C" void ScaLBL_D3Q7_AAeven_IonConcentration(double *dist, double *Den, int start, int finish, int Np){

	//cudaProfilerStart();
        /*
        DPCT1049:3: The work-group size passed to the SYCL kernel may exceed the
        limit. To get the device limit, query info::device::max_work_group_size.
        Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, NBLOCKS) *
                                  sycl::range<3>(1, 1, NTHREADS),
                              sycl::range<3>(1, 1, NTHREADS)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_AAeven_IonConcentration(
                        dist, Den, start, finish, Np, item_ct1);
            });

        /*
        DPCT1010:106: SYCL uses exceptions to report errors and does not use the
        error codes. The call was replaced with 0. You need to rewrite this
        code.
        */
        int err = 0;

        //cudaProfilerStop();
}

extern "C" void ScaLBL_D3Q7_AAodd_Ion(int *neighborList, double *dist, double *Den, double *FluxDiffusive, double *FluxAdvective, double *FluxElectrical, double *Velocity, double *ElectricField,  
                                      double Di, int zi, double rlx, double Vt, int start, int finish, int Np){
	//cudaProfilerStart();
        /*
        DPCT1049:4: The work-group size passed to the SYCL kernel may exceed the
        limit. To get the device limit, query info::device::max_work_group_size.
        Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, NBLOCKS) *
                                  sycl::range<3>(1, 1, NTHREADS),
                              sycl::range<3>(1, 1, NTHREADS)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_AAodd_Ion(
                        neighborList, dist, Den, FluxDiffusive, FluxAdvective,
                        FluxElectrical, Velocity, ElectricField, Di, zi, rlx,
                        Vt, start, finish, Np, item_ct1);
            });

        /*
        DPCT1010:108: SYCL uses exceptions to report errors and does not use the
        error codes. The call was replaced with 0. You need to rewrite this
        code.
        */
        int err = 0;

        //cudaProfilerStop();
}

extern "C" void ScaLBL_D3Q7_AAeven_Ion(double *dist, double *Den, double *FluxDiffusive, double *FluxAdvective, double *FluxElectrical, double *Velocity, double *ElectricField, 
                                       double Di, int zi, double rlx, double Vt, int start, int finish, int Np){
	//cudaProfilerStart();
        /*
        DPCT1049:5: The work-group size passed to the SYCL kernel may exceed the
        limit. To get the device limit, query info::device::max_work_group_size.
        Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, NBLOCKS) *
                                  sycl::range<3>(1, 1, NTHREADS),
                              sycl::range<3>(1, 1, NTHREADS)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_AAeven_Ion(
                        dist, Den, FluxDiffusive, FluxAdvective, FluxElectrical,
                        Velocity, ElectricField, Di, zi, rlx, Vt, start, finish,
                        Np, item_ct1);
            });

        /*
        DPCT1010:110: SYCL uses exceptions to report errors and does not use the
        error codes. The call was replaced with 0. You need to rewrite this
        code.
        */
        int err = 0;

        //cudaProfilerStop();
}

extern "C" void ScaLBL_D3Q7_Ion_Init(double *dist, double *Den, double DenInit, int Np){

	//cudaProfilerStart();
        /*
        DPCT1049:6: The work-group size passed to the SYCL kernel may exceed the
        limit. To get the device limit, query info::device::max_work_group_size.
        Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, NBLOCKS) *
                                  sycl::range<3>(1, 1, NTHREADS),
                              sycl::range<3>(1, 1, NTHREADS)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_Ion_Init(dist, Den, DenInit, Np, item_ct1);
            });

        /*
        DPCT1010:112: SYCL uses exceptions to report errors and does not use the
        error codes. The call was replaced with 0. You need to rewrite this
        code.
        */
        int err = 0;

        //cudaProfilerStop();
}

extern "C" void ScaLBL_D3Q7_Ion_Init_FromFile(double *dist, double *Den, int Np){

	//cudaProfilerStart();
        /*
        DPCT1049:7: The work-group size passed to the SYCL kernel may exceed the
        limit. To get the device limit, query info::device::max_work_group_size.
        Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, NBLOCKS) *
                                  sycl::range<3>(1, 1, NTHREADS),
                              sycl::range<3>(1, 1, NTHREADS)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_Ion_Init_FromFile(dist, Den, Np, item_ct1);
            });

        /*
        DPCT1010:114: SYCL uses exceptions to report errors and does not use the
        error codes. The call was replaced with 0. You need to rewrite this
        code.
        */
        int err = 0;

        //cudaProfilerStop();
}

extern "C" void ScaLBL_D3Q7_Ion_ChargeDensity(double *Den, double *ChargeDensity, double IonValence, int ion_component, int start, int finish, int Np){

	//cudaProfilerStart();
        /*
        DPCT1049:8: The work-group size passed to the SYCL kernel may exceed the
        limit. To get the device limit, query info::device::max_work_group_size.
        Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, NBLOCKS) *
                                  sycl::range<3>(1, 1, NTHREADS),
                              sycl::range<3>(1, 1, NTHREADS)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_Ion_ChargeDensity(
                        Den, ChargeDensity, IonValence, ion_component, start,
                        finish, Np, item_ct1);
            });

        /*
        DPCT1010:116: SYCL uses exceptions to report errors and does not use the
        error codes. The call was replaced with 0. You need to rewrite this
        code.
        */
        int err = 0;

        //cudaProfilerStop();
}

extern "C" void ScaLBL_D3Q7_Membrane_AssignLinkCoef(int *membrane, int *Map, double *Distance, double *Psi, double *coef,
		double Threshold, double MassFractionIn, double MassFractionOut, double ThresholdMassFractionIn, double ThresholdMassFractionOut,
		int memLinks, int Nx, int Ny, int Nz, int Np){

        /*
        DPCT1049:9: The work-group size passed to the SYCL kernel may exceed the
        limit. To get the device limit, query info::device::max_work_group_size.
        Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, NBLOCKS) *
                                  sycl::range<3>(1, 1, NTHREADS),
                              sycl::range<3>(1, 1, NTHREADS)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_Membrane_AssignLinkCoef(
                        membrane, Map, Distance, Psi, coef, Threshold,
                        MassFractionIn, MassFractionOut,
                        ThresholdMassFractionIn, ThresholdMassFractionOut,
                        memLinks, Nx, Ny, Nz, Np, item_ct1);
            });

        /*
        DPCT1010:118: SYCL uses exceptions to report errors and does not use the
        error codes. The call was replaced with 0. You need to rewrite this
        code.
        */
        int err = 0;
}

extern "C" void ScaLBL_D3Q7_Membrane_AssignLinkCoef_halo(
		const int Cqx, const int Cqy, int const Cqz, 
		int *Map, double *Distance, double *Psi, double Threshold, 
		double MassFractionIn, double MassFractionOut, double ThresholdMassFractionIn, double ThresholdMassFractionOut,
		int *d3q7_recvlist, int *d3q7_linkList, double *coef, int start, int nlinks, int count,
		const int N, const int Nx, const int Ny, const int Nz) {
	
    int GRID = count / NTHREADS + 1;

        /*
        DPCT1049:10: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, NTHREADS),
                              sycl::range<3>(1, 1, NTHREADS)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_Membrane_AssignLinkCoef_halo(
                        Cqx, Cqy, Cqz, Map, Distance, Psi, Threshold,
                        MassFractionIn, MassFractionOut,
                        ThresholdMassFractionIn, ThresholdMassFractionOut,
                        d3q7_recvlist, d3q7_linkList, coef, start, nlinks,
                        count, N, Nx, Ny, Nz, item_ct1);
            });

        /*
        DPCT1010:120: SYCL uses exceptions to report errors and does not use the
        error codes. The call was replaced with 0. You need to rewrite this
        code.
        */
        int err = 0;
}


extern "C" void ScaLBL_D3Q7_Membrane_Unpack(int q,  
		int *d3q7_recvlist, double *recvbuf, int count,
		double *dist, int N,  double *coef){
	
    int GRID = count / NTHREADS + 1;

        /*
        DPCT1049:11: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) *
                                  sycl::range<3>(1, 1, NTHREADS),
                              sycl::range<3>(1, 1, NTHREADS)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_Membrane_Unpack(q, d3q7_recvlist, recvbuf,
                                                    count, dist, N, coef,
                                                    item_ct1);
            });

        /*
        DPCT1010:122: SYCL uses exceptions to report errors and does not use the
        error codes. The call was replaced with 0. You need to rewrite this
        code.
        */
        int err = 0;
}

extern "C" void ScaLBL_D3Q7_Membrane_IonTransport(int *membrane, double *coef, 
		double *dist, double *Den, int memLinks, int Np){

        /*
        DPCT1049:12: The work-group size passed to the SYCL kernel may exceed
        the limit. To get the device limit, query
        info::device::max_work_group_size. Adjust the work-group size if needed.
        */
        dpct::get_default_queue().parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, NBLOCKS) *
                                  sycl::range<3>(1, 1, NTHREADS),
                              sycl::range<3>(1, 1, NTHREADS)),
            [=](sycl::nd_item<3> item_ct1) {
                    dvc_ScaLBL_D3Q7_Membrane_IonTransport(
                        membrane, coef, dist, Den, memLinks, Np, item_ct1);
            });

        /*
        DPCT1010:124: SYCL uses exceptions to report errors and does not use the
        error codes. The call was replaced with 0. You need to rewrite this
        code.
        */
        int err = 0;
}


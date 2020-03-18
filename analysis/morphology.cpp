#include <analysis/morphology.h>
// Implementation of morphological opening routine

inline void PackID(int *list, int count, signed char *sendbuf, signed char *ID){
	// Fill in the phase ID values from neighboring processors
	// This packs up the values that need to be sent from one processor to another
	int idx,n;

	for (idx=0; idx<count; idx++){
		n = list[idx];
		sendbuf[idx] = ID[n];
	}
}
//***************************************************************************************

inline void UnpackID(int *list, int count, signed char *recvbuf, signed char *ID){
	// Fill in the phase ID values from neighboring processors
	// This unpacks the values once they have been recieved from neighbors
	int idx,n;

	for (idx=0; idx<count; idx++){
		n = list[idx];
		ID[n] = recvbuf[idx];
	}
}

//***************************************************************************************
double MorphOpen(DoubleArray &SignDist, signed char *id, std::shared_ptr<Domain> Dm, double VoidFraction, signed char ErodeLabel, signed char NewLabel){
	// SignDist is the distance to the object that you want to constaing the morphological opening
	// VoidFraction is the the empty space where the object inst
	// id is a labeled map
	// Dm contains information about the domain structure
	
	int nx = Dm->Nx;
	int ny = Dm->Ny;
	int nz = Dm->Nz;
	int nprocx = Dm->nprocx();
	int nprocy = Dm->nprocy();
	int nprocz = Dm->nprocz();
	int rank = Dm->rank();

	int n;
	double final_void_fraction;
	double count,countGlobal,totalGlobal;
	count = 0.f;
	double maxdist=-200.f;
	double maxdistGlobal;
	for (int k=1; k<nz-1; k++){
		for (int j=1; j<ny-1; j++){
			for (int i=1; i<nx-1; i++){
				n = k*nx*ny+j*nx+i;
				// extract maximum distance for critical radius
				if ( SignDist(i,j,k) > maxdist) maxdist=SignDist(i,j,k);
				if ( id[n] == ErodeLabel){
					count += 1.0;
					//id[n]  = 2;
				}
			}
		}
	}
	Dm->Comm.barrier();
	
	// total Global is the number of nodes in the pore-space
	totalGlobal = Dm->Comm.sumReduce( count );
	maxdistGlobal = Dm->Comm.sumReduce( maxdist );
	double volume=double(nprocx*nprocy*nprocz)*double(nx-2)*double(ny-2)*double(nz-2);
	double volume_fraction=totalGlobal/volume;
	if (rank==0) printf("Volume fraction for morphological opening: %f \n",volume_fraction);
	if (rank==0) printf("Maximum pore size: %f \n",maxdistGlobal);
	final_void_fraction = volume_fraction; //initialize

	// Communication buffers
	signed char *sendID_x, *sendID_y, *sendID_z, *sendID_X, *sendID_Y, *sendID_Z;
	signed char *sendID_xy, *sendID_yz, *sendID_xz, *sendID_Xy, *sendID_Yz, *sendID_xZ;
	signed char *sendID_xY, *sendID_yZ, *sendID_Xz, *sendID_XY, *sendID_YZ, *sendID_XZ;
	signed char *recvID_x, *recvID_y, *recvID_z, *recvID_X, *recvID_Y, *recvID_Z;
	signed char *recvID_xy, *recvID_yz, *recvID_xz, *recvID_Xy, *recvID_Yz, *recvID_xZ;
	signed char *recvID_xY, *recvID_yZ, *recvID_Xz, *recvID_XY, *recvID_YZ, *recvID_XZ;
	// send buffers
	sendID_x = new signed char [Dm->sendCount_x];
	sendID_y = new signed char [Dm->sendCount_y];
	sendID_z = new signed char [Dm->sendCount_z];
	sendID_X = new signed char [Dm->sendCount_X];
	sendID_Y = new signed char [Dm->sendCount_Y];
	sendID_Z = new signed char [Dm->sendCount_Z];
	sendID_xy = new signed char [Dm->sendCount_xy];
	sendID_yz = new signed char [Dm->sendCount_yz];
	sendID_xz = new signed char [Dm->sendCount_xz];
	sendID_Xy = new signed char [Dm->sendCount_Xy];
	sendID_Yz = new signed char [Dm->sendCount_Yz];
	sendID_xZ = new signed char [Dm->sendCount_xZ];
	sendID_xY = new signed char [Dm->sendCount_xY];
	sendID_yZ = new signed char [Dm->sendCount_yZ];
	sendID_Xz = new signed char [Dm->sendCount_Xz];
	sendID_XY = new signed char [Dm->sendCount_XY];
	sendID_YZ = new signed char [Dm->sendCount_YZ];
	sendID_XZ = new signed char [Dm->sendCount_XZ];
	//......................................................................................
	// recv buffers
	recvID_x = new signed char [Dm->recvCount_x];
	recvID_y = new signed char [Dm->recvCount_y];
	recvID_z = new signed char [Dm->recvCount_z];
	recvID_X = new signed char [Dm->recvCount_X];
	recvID_Y = new signed char [Dm->recvCount_Y];
	recvID_Z = new signed char [Dm->recvCount_Z];
	recvID_xy = new signed char [Dm->recvCount_xy];
	recvID_yz = new signed char [Dm->recvCount_yz];
	recvID_xz = new signed char [Dm->recvCount_xz];
	recvID_Xy = new signed char [Dm->recvCount_Xy];
	recvID_xZ = new signed char [Dm->recvCount_xZ];
	recvID_xY = new signed char [Dm->recvCount_xY];
	recvID_yZ = new signed char [Dm->recvCount_yZ];
	recvID_Yz = new signed char [Dm->recvCount_Yz];
	recvID_Xz = new signed char [Dm->recvCount_Xz];
	recvID_XY = new signed char [Dm->recvCount_XY];
	recvID_YZ = new signed char [Dm->recvCount_YZ];
	recvID_XZ = new signed char [Dm->recvCount_XZ];
	//......................................................................................
	int sendtag,recvtag;
	sendtag = recvtag = 7;

	int ii,jj,kk;
	int Nx = nx;
	int Ny = ny;
	int Nz = nz;

	double void_fraction_old=1.0;
	double void_fraction_new=1.0; 
	double void_fraction_diff_old = 1.0;
	double void_fraction_diff_new = 1.0;

	// Increase the critical radius until the target saturation is met
	double deltaR=0.05; // amount to change the radius in voxel units
	double Rcrit_old=0.0;

	int imin,jmin,kmin,imax,jmax,kmax;

	if (ErodeLabel == 1){
		VoidFraction = 1.0 - VoidFraction;
	}

	double Rcrit_new = maxdistGlobal;

	while (void_fraction_new > VoidFraction)
	{
		void_fraction_diff_old = void_fraction_diff_new;
		void_fraction_old = void_fraction_new;
		Rcrit_old = Rcrit_new;
		Rcrit_new -= deltaR*Rcrit_old;
		int Window=round(Rcrit_new);
		if (Window == 0) Window = 1; // If Window = 0 at the begining, after the following process will have sw=1.0
		// and sw<Sw will be immediately broken
		double LocalNumber=0.f;
		for(int k=0; k<Nz; k++){
			for(int j=0; j<Ny; j++){
				for(int i=0; i<Nx; i++){
					n = k*nx*ny + j*nx+i;
					if (SignDist(i,j,k) > Rcrit_new){
						// loop over the window and update
						imin=max(1,i-Window);
						jmin=max(1,j-Window);
						kmin=max(1,k-Window);
						imax=min(Nx-1,i+Window);
						jmax=min(Ny-1,j+Window);
						kmax=min(Nz-1,k+Window);
						for (kk=kmin; kk<kmax; kk++){
							for (jj=jmin; jj<jmax; jj++){
								for (ii=imin; ii<imax; ii++){
									int nn = kk*nx*ny+jj*nx+ii;
									double dsq = double((ii-i)*(ii-i)+(jj-j)*(jj-j)+(kk-k)*(kk-k));
									if (id[nn] == ErodeLabel && dsq <= Rcrit_new*Rcrit_new){
										LocalNumber+=1.0;
										id[nn]=NewLabel;
									}
								}
							}
						}

					}
					// move on
				}
			}
		}
		// Pack and send the updated ID values
		PackID(Dm->sendList_x, Dm->sendCount_x ,sendID_x, id);
		PackID(Dm->sendList_X, Dm->sendCount_X ,sendID_X, id);
		PackID(Dm->sendList_y, Dm->sendCount_y ,sendID_y, id);
		PackID(Dm->sendList_Y, Dm->sendCount_Y ,sendID_Y, id);
		PackID(Dm->sendList_z, Dm->sendCount_z ,sendID_z, id);
		PackID(Dm->sendList_Z, Dm->sendCount_Z ,sendID_Z, id);
		PackID(Dm->sendList_xy, Dm->sendCount_xy ,sendID_xy, id);
		PackID(Dm->sendList_Xy, Dm->sendCount_Xy ,sendID_Xy, id);
		PackID(Dm->sendList_xY, Dm->sendCount_xY ,sendID_xY, id);
		PackID(Dm->sendList_XY, Dm->sendCount_XY ,sendID_XY, id);
		PackID(Dm->sendList_xz, Dm->sendCount_xz ,sendID_xz, id);
		PackID(Dm->sendList_Xz, Dm->sendCount_Xz ,sendID_Xz, id);
		PackID(Dm->sendList_xZ, Dm->sendCount_xZ ,sendID_xZ, id);
		PackID(Dm->sendList_XZ, Dm->sendCount_XZ ,sendID_XZ, id);
		PackID(Dm->sendList_yz, Dm->sendCount_yz ,sendID_yz, id);
		PackID(Dm->sendList_Yz, Dm->sendCount_Yz ,sendID_Yz, id);
		PackID(Dm->sendList_yZ, Dm->sendCount_yZ ,sendID_yZ, id);
		PackID(Dm->sendList_YZ, Dm->sendCount_YZ ,sendID_YZ, id);
		//......................................................................................
		MPI_Sendrecv(sendID_x,Dm->sendCount_x,MPI_CHAR,Dm->rank_x(),sendtag,
				recvID_X,Dm->recvCount_X,MPI_CHAR,Dm->rank_X(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_X,Dm->sendCount_X,MPI_CHAR,Dm->rank_X(),sendtag,
				recvID_x,Dm->recvCount_x,MPI_CHAR,Dm->rank_x(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_y,Dm->sendCount_y,MPI_CHAR,Dm->rank_y(),sendtag,
				recvID_Y,Dm->recvCount_Y,MPI_CHAR,Dm->rank_Y(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_Y,Dm->sendCount_Y,MPI_CHAR,Dm->rank_Y(),sendtag,
				recvID_y,Dm->recvCount_y,MPI_CHAR,Dm->rank_y(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_z,Dm->sendCount_z,MPI_CHAR,Dm->rank_z(),sendtag,
				recvID_Z,Dm->recvCount_Z,MPI_CHAR,Dm->rank_Z(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_Z,Dm->sendCount_Z,MPI_CHAR,Dm->rank_Z(),sendtag,
				recvID_z,Dm->recvCount_z,MPI_CHAR,Dm->rank_z(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_xy,Dm->sendCount_xy,MPI_CHAR,Dm->rank_xy(),sendtag,
				recvID_XY,Dm->recvCount_XY,MPI_CHAR,Dm->rank_XY(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_XY,Dm->sendCount_XY,MPI_CHAR,Dm->rank_XY(),sendtag,
				recvID_xy,Dm->recvCount_xy,MPI_CHAR,Dm->rank_xy(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_Xy,Dm->sendCount_Xy,MPI_CHAR,Dm->rank_Xy(),sendtag,
				recvID_xY,Dm->recvCount_xY,MPI_CHAR,Dm->rank_xY(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_xY,Dm->sendCount_xY,MPI_CHAR,Dm->rank_xY(),sendtag,
				recvID_Xy,Dm->recvCount_Xy,MPI_CHAR,Dm->rank_Xy(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_xz,Dm->sendCount_xz,MPI_CHAR,Dm->rank_xz(),sendtag,
				recvID_XZ,Dm->recvCount_XZ,MPI_CHAR,Dm->rank_XZ(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_XZ,Dm->sendCount_XZ,MPI_CHAR,Dm->rank_XZ(),sendtag,
				recvID_xz,Dm->recvCount_xz,MPI_CHAR,Dm->rank_xz(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_Xz,Dm->sendCount_Xz,MPI_CHAR,Dm->rank_Xz(),sendtag,
				recvID_xZ,Dm->recvCount_xZ,MPI_CHAR,Dm->rank_xZ(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_xZ,Dm->sendCount_xZ,MPI_CHAR,Dm->rank_xZ(),sendtag,
				recvID_Xz,Dm->recvCount_Xz,MPI_CHAR,Dm->rank_Xz(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_yz,Dm->sendCount_yz,MPI_CHAR,Dm->rank_yz(),sendtag,
				recvID_YZ,Dm->recvCount_YZ,MPI_CHAR,Dm->rank_YZ(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_YZ,Dm->sendCount_YZ,MPI_CHAR,Dm->rank_YZ(),sendtag,
				recvID_yz,Dm->recvCount_yz,MPI_CHAR,Dm->rank_yz(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_Yz,Dm->sendCount_Yz,MPI_CHAR,Dm->rank_Yz(),sendtag,
				recvID_yZ,Dm->recvCount_yZ,MPI_CHAR,Dm->rank_yZ(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_yZ,Dm->sendCount_yZ,MPI_CHAR,Dm->rank_yZ(),sendtag,
				recvID_Yz,Dm->recvCount_Yz,MPI_CHAR,Dm->rank_Yz(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		//......................................................................................
		UnpackID(Dm->recvList_x, Dm->recvCount_x ,recvID_x, id);
		UnpackID(Dm->recvList_X, Dm->recvCount_X ,recvID_X, id);
		UnpackID(Dm->recvList_y, Dm->recvCount_y ,recvID_y, id);
		UnpackID(Dm->recvList_Y, Dm->recvCount_Y ,recvID_Y, id);
		UnpackID(Dm->recvList_z, Dm->recvCount_z ,recvID_z, id);
		UnpackID(Dm->recvList_Z, Dm->recvCount_Z ,recvID_Z, id);
		UnpackID(Dm->recvList_xy, Dm->recvCount_xy ,recvID_xy, id);
		UnpackID(Dm->recvList_Xy, Dm->recvCount_Xy ,recvID_Xy, id);
		UnpackID(Dm->recvList_xY, Dm->recvCount_xY ,recvID_xY, id);
		UnpackID(Dm->recvList_XY, Dm->recvCount_XY ,recvID_XY, id);
		UnpackID(Dm->recvList_xz, Dm->recvCount_xz ,recvID_xz, id);
		UnpackID(Dm->recvList_Xz, Dm->recvCount_Xz ,recvID_Xz, id);
		UnpackID(Dm->recvList_xZ, Dm->recvCount_xZ ,recvID_xZ, id);
		UnpackID(Dm->recvList_XZ, Dm->recvCount_XZ ,recvID_XZ, id);
		UnpackID(Dm->recvList_yz, Dm->recvCount_yz ,recvID_yz, id);
		UnpackID(Dm->recvList_Yz, Dm->recvCount_Yz ,recvID_Yz, id);
		UnpackID(Dm->recvList_yZ, Dm->recvCount_yZ ,recvID_yZ, id);
		UnpackID(Dm->recvList_YZ, Dm->recvCount_YZ ,recvID_YZ, id);
		//......................................................................................

		//double GlobalNumber = Dm->Comm.sumReduce( LocalNumber );

		count = 0.f;
		for (int k=1; k<Nz-1; k++){
			for (int j=1; j<Ny-1; j++){
				for (int i=1; i<Nx-1; i++){
					n=k*Nx*Ny+j*Nx+i;
					if (id[n] == ErodeLabel){
						count+=1.0;
					}
				}
			}
		}
		countGlobal = Dm->Comm.sumReduce( count );
		void_fraction_new = countGlobal/totalGlobal;
		void_fraction_diff_new = abs(void_fraction_new-VoidFraction);
	/*	if (rank==0){
			printf("     %f ",void_fraction_new);
			printf("     %f\n",Rcrit_new);
		} */
	}

	if (void_fraction_diff_new<void_fraction_diff_old){
		final_void_fraction=void_fraction_new;
		if (rank==0){
			printf("Final void fraction =%f\n",void_fraction_new);
			printf("Final critical radius=%f\n",Rcrit_new);
		}
	}
	else{
		final_void_fraction=void_fraction_old;
		if (rank==0){
			printf("Final void fraction=%f\n",void_fraction_old);
			printf("Final critical radius=%f\n",Rcrit_old);
		}
	}
	return final_void_fraction;
}
/*
double morph_open()
{

	fillHalo<char> fillChar(Dm->Comm,Dm->rank_info,{Nx-2,Ny-2,Nz-2},{1,1,1},0,1);


	MPI_Allreduce(&LocalNumber,&GlobalNumber,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);

	count = 0.f;
	for (int k=1; k<Nz-1; k++){
		for (int j=1; j<Ny-1; j++){
			for (int i=1; i<Nx-1; i++){
				n=k*Nx*Ny+j*Nx+i;
				if (id[n] == 2){
					count+=1.0;
				}
			}
		}
	}
	MPI_Allreduce(&count,&countGlobal,1,MPI_DOUBLE,MPI_SUM,Dm->Comm);
	return countGlobal;
}
*/

//***************************************************************************************
double MorphDrain(DoubleArray &SignDist, signed char *id, std::shared_ptr<Domain> Dm, double VoidFraction){
	// SignDist is the distance to the object that you want to constaing the morphological opening
	// VoidFraction is the the empty space where the object inst
	// id is a labeled map
	// Dm contains information about the domain structure
	
	int nx = Dm->Nx;
	int ny = Dm->Ny;
	int nz = Dm->Nz;
	int nprocx = Dm->nprocx();
	int nprocy = Dm->nprocy();
	int nprocz = Dm->nprocz();
	int rank = Dm->rank();
	
	DoubleArray phase(nx,ny,nz);
	IntArray phase_label(nx,ny,nz);

	int n;
	double final_void_fraction;
	double count,countGlobal,totalGlobal;
	count = 0.f;
	double maxdist=-200.f;
	double maxdistGlobal;
	for (int k=1; k<nz-1; k++){
		for (int j=1; j<ny-1; j++){
			for (int i=1; i<nx-1; i++){
				n = k*nx*ny+j*nx+i;
				// extract maximum distance for critical radius
				if ( SignDist(i,j,k) > maxdist) maxdist=SignDist(i,j,k);
				if ( SignDist(i,j,k) > 0.0 ){
					count += 1.0;
					id[n]  = 2;
				}
			}
		}
	}

	Dm->Comm.barrier();
	
	// total Global is the number of nodes in the pore-space
	totalGlobal = Dm->Comm.sumReduce( count );
	maxdistGlobal = Dm->Comm.sumReduce( maxdist );
	double volume=double(nprocx*nprocy*nprocz)*double(nx-2)*double(ny-2)*double(nz-2);
	double volume_fraction=totalGlobal/volume;
	if (rank==0) printf("Volume fraction for morphological opening: %f \n",volume_fraction);
	if (rank==0) printf("Maximum pore size: %f \n",maxdistGlobal);

	// Communication buffers
	signed char *sendID_x, *sendID_y, *sendID_z, *sendID_X, *sendID_Y, *sendID_Z;
	signed char *sendID_xy, *sendID_yz, *sendID_xz, *sendID_Xy, *sendID_Yz, *sendID_xZ;
	signed char *sendID_xY, *sendID_yZ, *sendID_Xz, *sendID_XY, *sendID_YZ, *sendID_XZ;
	signed char *recvID_x, *recvID_y, *recvID_z, *recvID_X, *recvID_Y, *recvID_Z;
	signed char *recvID_xy, *recvID_yz, *recvID_xz, *recvID_Xy, *recvID_Yz, *recvID_xZ;
	signed char *recvID_xY, *recvID_yZ, *recvID_Xz, *recvID_XY, *recvID_YZ, *recvID_XZ;
	// send buffers
	sendID_x = new signed char [Dm->sendCount_x];
	sendID_y = new signed char [Dm->sendCount_y];
	sendID_z = new signed char [Dm->sendCount_z];
	sendID_X = new signed char [Dm->sendCount_X];
	sendID_Y = new signed char [Dm->sendCount_Y];
	sendID_Z = new signed char [Dm->sendCount_Z];
	sendID_xy = new signed char [Dm->sendCount_xy];
	sendID_yz = new signed char [Dm->sendCount_yz];
	sendID_xz = new signed char [Dm->sendCount_xz];
	sendID_Xy = new signed char [Dm->sendCount_Xy];
	sendID_Yz = new signed char [Dm->sendCount_Yz];
	sendID_xZ = new signed char [Dm->sendCount_xZ];
	sendID_xY = new signed char [Dm->sendCount_xY];
	sendID_yZ = new signed char [Dm->sendCount_yZ];
	sendID_Xz = new signed char [Dm->sendCount_Xz];
	sendID_XY = new signed char [Dm->sendCount_XY];
	sendID_YZ = new signed char [Dm->sendCount_YZ];
	sendID_XZ = new signed char [Dm->sendCount_XZ];
	//......................................................................................
	// recv buffers
	recvID_x = new signed char [Dm->recvCount_x];
	recvID_y = new signed char [Dm->recvCount_y];
	recvID_z = new signed char [Dm->recvCount_z];
	recvID_X = new signed char [Dm->recvCount_X];
	recvID_Y = new signed char [Dm->recvCount_Y];
	recvID_Z = new signed char [Dm->recvCount_Z];
	recvID_xy = new signed char [Dm->recvCount_xy];
	recvID_yz = new signed char [Dm->recvCount_yz];
	recvID_xz = new signed char [Dm->recvCount_xz];
	recvID_Xy = new signed char [Dm->recvCount_Xy];
	recvID_xZ = new signed char [Dm->recvCount_xZ];
	recvID_xY = new signed char [Dm->recvCount_xY];
	recvID_yZ = new signed char [Dm->recvCount_yZ];
	recvID_Yz = new signed char [Dm->recvCount_Yz];
	recvID_Xz = new signed char [Dm->recvCount_Xz];
	recvID_XY = new signed char [Dm->recvCount_XY];
	recvID_YZ = new signed char [Dm->recvCount_YZ];
	recvID_XZ = new signed char [Dm->recvCount_XZ];
	//......................................................................................
	int sendtag,recvtag;
	sendtag = recvtag = 7;

	int ii,jj,kk;
	int Nx = nx;
	int Ny = ny;
	int Nz = nz;

	double void_fraction_old=1.0;
	double void_fraction_new=1.0; 
	double void_fraction_diff_old = 1.0;
	double void_fraction_diff_new = 1.0;

	// Increase the critical radius until the target saturation is met
	double deltaR=0.05; // amount to change the radius in voxel units
	double Rcrit_old;

	int imin,jmin,kmin,imax,jmax,kmax;

	double Rcrit_new = maxdistGlobal;
	//if (argc>2){
	//	Rcrit_new = strtod(argv[2],NULL);
	//	if (rank==0) printf("Max. distance =%f, Initial critical radius = %f \n",maxdistGlobal,Rcrit_new);
	//}
	Dm->Comm.barrier();

	
	FILE *DRAIN = fopen("morphdrain.csv","w");
	fprintf(DRAIN,"sw radius\n");				

	while (void_fraction_new > VoidFraction && Rcrit_new > 0.5)
	{
		void_fraction_diff_old = void_fraction_diff_new;
		void_fraction_old = void_fraction_new;
		Rcrit_old = Rcrit_new;
		Rcrit_new -= deltaR*Rcrit_old;
		int Window=round(Rcrit_new);
		if (Window == 0) Window = 1; // If Window = 0 at the begining, after the following process will have sw=1.0
		// and sw<Sw will be immediately broken
		double LocalNumber=0.f;
		for(int k=1; k<Nz-1; k++){
			for(int j=1; j<Ny-1; j++){
				for(int i=1; i<Nx-1; i++){
					n = k*nx*ny + j*nx+i;
					if (SignDist(i,j,k) > Rcrit_new){
						// loop over the window and update
						imin=max(1,i-Window);
						jmin=max(1,j-Window);
						kmin=max(1,k-Window);
						imax=min(Nx-1,i+Window);
						jmax=min(Ny-1,j+Window);
						kmax=min(Nz-1,k+Window);
						for (kk=kmin; kk<kmax; kk++){
							for (jj=jmin; jj<jmax; jj++){
								for (ii=imin; ii<imax; ii++){
									int nn = kk*nx*ny+jj*nx+ii;
									double dsq = double((ii-i)*(ii-i)+(jj-j)*(jj-j)+(kk-k)*(kk-k));
									if (id[nn] == 2 && dsq <= (Rcrit_new+1)*(Rcrit_new+1)){
										LocalNumber+=1.0;
										id[nn]=1;
									}
								}
							}
						}

					}
					// move on
				}
			}
		}
		// Pack and send the updated ID values
		PackID(Dm->sendList_x, Dm->sendCount_x ,sendID_x, id);
		PackID(Dm->sendList_X, Dm->sendCount_X ,sendID_X, id);
		PackID(Dm->sendList_y, Dm->sendCount_y ,sendID_y, id);
		PackID(Dm->sendList_Y, Dm->sendCount_Y ,sendID_Y, id);
		PackID(Dm->sendList_z, Dm->sendCount_z ,sendID_z, id);
		PackID(Dm->sendList_Z, Dm->sendCount_Z ,sendID_Z, id);
		PackID(Dm->sendList_xy, Dm->sendCount_xy ,sendID_xy, id);
		PackID(Dm->sendList_Xy, Dm->sendCount_Xy ,sendID_Xy, id);
		PackID(Dm->sendList_xY, Dm->sendCount_xY ,sendID_xY, id);
		PackID(Dm->sendList_XY, Dm->sendCount_XY ,sendID_XY, id);
		PackID(Dm->sendList_xz, Dm->sendCount_xz ,sendID_xz, id);
		PackID(Dm->sendList_Xz, Dm->sendCount_Xz ,sendID_Xz, id);
		PackID(Dm->sendList_xZ, Dm->sendCount_xZ ,sendID_xZ, id);
		PackID(Dm->sendList_XZ, Dm->sendCount_XZ ,sendID_XZ, id);
		PackID(Dm->sendList_yz, Dm->sendCount_yz ,sendID_yz, id);
		PackID(Dm->sendList_Yz, Dm->sendCount_Yz ,sendID_Yz, id);
		PackID(Dm->sendList_yZ, Dm->sendCount_yZ ,sendID_yZ, id);
		PackID(Dm->sendList_YZ, Dm->sendCount_YZ ,sendID_YZ, id);
		//......................................................................................
		MPI_Sendrecv(sendID_x,Dm->sendCount_x,MPI_CHAR,Dm->rank_x(),sendtag,
				recvID_X,Dm->recvCount_X,MPI_CHAR,Dm->rank_X(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_X,Dm->sendCount_X,MPI_CHAR,Dm->rank_X(),sendtag,
				recvID_x,Dm->recvCount_x,MPI_CHAR,Dm->rank_x(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_y,Dm->sendCount_y,MPI_CHAR,Dm->rank_y(),sendtag,
				recvID_Y,Dm->recvCount_Y,MPI_CHAR,Dm->rank_Y(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_Y,Dm->sendCount_Y,MPI_CHAR,Dm->rank_Y(),sendtag,
				recvID_y,Dm->recvCount_y,MPI_CHAR,Dm->rank_y(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_z,Dm->sendCount_z,MPI_CHAR,Dm->rank_z(),sendtag,
				recvID_Z,Dm->recvCount_Z,MPI_CHAR,Dm->rank_Z(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_Z,Dm->sendCount_Z,MPI_CHAR,Dm->rank_Z(),sendtag,
				recvID_z,Dm->recvCount_z,MPI_CHAR,Dm->rank_z(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_xy,Dm->sendCount_xy,MPI_CHAR,Dm->rank_xy(),sendtag,
				recvID_XY,Dm->recvCount_XY,MPI_CHAR,Dm->rank_XY(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_XY,Dm->sendCount_XY,MPI_CHAR,Dm->rank_XY(),sendtag,
				recvID_xy,Dm->recvCount_xy,MPI_CHAR,Dm->rank_xy(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_Xy,Dm->sendCount_Xy,MPI_CHAR,Dm->rank_Xy(),sendtag,
				recvID_xY,Dm->recvCount_xY,MPI_CHAR,Dm->rank_xY(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_xY,Dm->sendCount_xY,MPI_CHAR,Dm->rank_xY(),sendtag,
				recvID_Xy,Dm->recvCount_Xy,MPI_CHAR,Dm->rank_Xy(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_xz,Dm->sendCount_xz,MPI_CHAR,Dm->rank_xz(),sendtag,
				recvID_XZ,Dm->recvCount_XZ,MPI_CHAR,Dm->rank_XZ(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_XZ,Dm->sendCount_XZ,MPI_CHAR,Dm->rank_XZ(),sendtag,
				recvID_xz,Dm->recvCount_xz,MPI_CHAR,Dm->rank_xz(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_Xz,Dm->sendCount_Xz,MPI_CHAR,Dm->rank_Xz(),sendtag,
				recvID_xZ,Dm->recvCount_xZ,MPI_CHAR,Dm->rank_xZ(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_xZ,Dm->sendCount_xZ,MPI_CHAR,Dm->rank_xZ(),sendtag,
				recvID_Xz,Dm->recvCount_Xz,MPI_CHAR,Dm->rank_Xz(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_yz,Dm->sendCount_yz,MPI_CHAR,Dm->rank_yz(),sendtag,
				recvID_YZ,Dm->recvCount_YZ,MPI_CHAR,Dm->rank_YZ(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_YZ,Dm->sendCount_YZ,MPI_CHAR,Dm->rank_YZ(),sendtag,
				recvID_yz,Dm->recvCount_yz,MPI_CHAR,Dm->rank_yz(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_Yz,Dm->sendCount_Yz,MPI_CHAR,Dm->rank_Yz(),sendtag,
				recvID_yZ,Dm->recvCount_yZ,MPI_CHAR,Dm->rank_yZ(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		MPI_Sendrecv(sendID_yZ,Dm->sendCount_yZ,MPI_CHAR,Dm->rank_yZ(),sendtag,
				recvID_Yz,Dm->recvCount_Yz,MPI_CHAR,Dm->rank_Yz(),recvtag,Dm->Comm.getCommunicator(),MPI_STATUS_IGNORE);
		//......................................................................................
		UnpackID(Dm->recvList_x, Dm->recvCount_x ,recvID_x, id);
		UnpackID(Dm->recvList_X, Dm->recvCount_X ,recvID_X, id);
		UnpackID(Dm->recvList_y, Dm->recvCount_y ,recvID_y, id);
		UnpackID(Dm->recvList_Y, Dm->recvCount_Y ,recvID_Y, id);
		UnpackID(Dm->recvList_z, Dm->recvCount_z ,recvID_z, id);
		UnpackID(Dm->recvList_Z, Dm->recvCount_Z ,recvID_Z, id);
		UnpackID(Dm->recvList_xy, Dm->recvCount_xy ,recvID_xy, id);
		UnpackID(Dm->recvList_Xy, Dm->recvCount_Xy ,recvID_Xy, id);
		UnpackID(Dm->recvList_xY, Dm->recvCount_xY ,recvID_xY, id);
		UnpackID(Dm->recvList_XY, Dm->recvCount_XY ,recvID_XY, id);
		UnpackID(Dm->recvList_xz, Dm->recvCount_xz ,recvID_xz, id);
		UnpackID(Dm->recvList_Xz, Dm->recvCount_Xz ,recvID_Xz, id);
		UnpackID(Dm->recvList_xZ, Dm->recvCount_xZ ,recvID_xZ, id);
		UnpackID(Dm->recvList_XZ, Dm->recvCount_XZ ,recvID_XZ, id);
		UnpackID(Dm->recvList_yz, Dm->recvCount_yz ,recvID_yz, id);
		UnpackID(Dm->recvList_Yz, Dm->recvCount_Yz ,recvID_Yz, id);
		UnpackID(Dm->recvList_yZ, Dm->recvCount_yZ ,recvID_yZ, id);
		UnpackID(Dm->recvList_YZ, Dm->recvCount_YZ ,recvID_YZ, id);
		//......................................................................................
		// double GlobalNumber = Dm->Comm.sumReduce( LocalNumber );
		
		for (int k=0; k<nz; k++){
			for (int j=0; j<ny; j++){
				for (int i=0; i<nx; i++){
					n=k*nx*ny+j*nx+i;
					if (id[n] == 1){
						phase(i,j,k) = 1.0;
					}
					else
						phase(i,j,k) = -1.0;
				}
			}
		}
		
		// Extract only the connected part of NWP
		BlobIDstruct new_index;
		double vF=0.0; double vS=0.0;
		ComputeGlobalBlobIDs(nx-2,ny-2,nz-2,Dm->rank_info,phase,SignDist,vF,vS,phase_label,Dm->Comm);
		Dm->Comm.barrier();
		
		for (int k=0; k<nz; k++){
			for (int j=0; j<ny; j++){
				for (int i=0; i<nx; i++){
					n=k*nx*ny+j*nx+i;
					if (id[n] == 1 && phase_label(i,j,k) > 1){
						id[n] = 2;
					}
				}
			}
		}
		
		
		/*
		* Extract only the connected part of NWP
		for (int k=1; k<nz-1; k++){
			for (int j=1; j<ny-1; j++){
				for (int i=1; i<nx-1; i++){
					n=k*nx*ny+j*nx+i;
					if (id[n] == 2){
						phase(i,j,k) = 1.0;
					}
					else if (id[n] == 1){
						// nwp
						phase(i,j,k) = -1.0;
					}
					else{
						// treat solid as WP since films can connect 
						phase(i,j,k) = 1.0;
					}
				}
			}
		}
		
		ComputeGlobalBlobIDs(nx-2,ny-2,nz-2,Dm->rank_info,phase,SignDist,vF,vS,phase_label,Dm->Comm);
		MPI_Barrier(Dm->Comm);
		
		for (int k=1; k<nz-1; k++){
			for (int j=1; j<ny-1; j++){
				for (int i=1; i<nx-1; i++){
					n=k*nx*ny+j*nx+i;
					if (id[n] == 2 && phase_label(i,j,k) > 1){
						id[n] = 20;
					}
				}
			}
		}
		// done
		*/

		count = 0.f;
		for (int k=1; k<nz-1; k++){
			for (int j=1; j<ny-1; j++){
				for (int i=1; i<nx-1; i++){
					n=k*nx*ny+j*nx+i;
					if (id[n] > 1){
						count+=1.0;
					}
				}
			}
		}
		countGlobal = Dm->Comm.sumReduce( count );
		void_fraction_new = countGlobal/totalGlobal;
		void_fraction_diff_new = abs(void_fraction_new-VoidFraction);
		if (rank==0){
			fprintf(DRAIN,"%f ",void_fraction_new);
			fprintf(DRAIN,"%f\n",Rcrit_new);
			printf("     %f ",void_fraction_new);
			printf("     %f\n",Rcrit_new);
		}
	}

	if (void_fraction_diff_new<void_fraction_diff_old){
		final_void_fraction=void_fraction_new;
		if (rank==0){
			printf("Final void fraction =%f\n",void_fraction_new);
			printf("Final critical radius=%f\n",Rcrit_new);
		}
	}
	else{
		final_void_fraction=void_fraction_old;
		if (rank==0){
			printf("Final void fraction=%f\n",void_fraction_old);
			printf("Final critical radius=%f\n",Rcrit_old);
		}
	}
	// label all WP components as 2
	for (int k=1; k<nz-1; k++){
		for (int j=1; j<ny-1; j++){
			for (int i=1; i<nx-1; i++){
				n=k*nx*ny+j*nx+i;
				if (id[n] > 1){
					id[n] = 2;
				}
			}
		}
	}
	
	return final_void_fraction;
}

double MorphGrow(DoubleArray &BoundaryDist, DoubleArray &Dist, Array<char> &id, std::shared_ptr<Domain> Dm, double TargetGrowth)
{
	int Nx = Dm->Nx;
	int Ny = Dm->Ny;
	int Nz = Dm->Nz;
	int rank = Dm->rank();
	
	double count=0.0;
	for (int k=1; k<Nz-1; k++){
		for (int j=1; j<Ny-1; j++){
			for (int i=1; i<Nx-1; i++){
				if (Dist(i,j,k) < 0.0){
					count+=1.0;
				}
			}
		}
	}
	double count_original = Dm->Comm.sumReduce( count);

	// Estimate morph_delta
	double morph_delta = 0.0;
	if (TargetGrowth > 0.0) morph_delta = 0.1;
	else					morph_delta = -0.1;
	double morph_delta_previous = 0.0;
	double GrowthEstimate = 0.0;
	double GrowthPrevious = 0.0;
	int COUNT_FOR_LOOP = 0;
	double ERROR = 100.0;
	if (rank == 0) printf("Estimate delta for growth=%f \n",TargetGrowth);
	while ( ERROR > 0.01 && COUNT_FOR_LOOP < 10 ){
		COUNT_FOR_LOOP++;
		count = 0.0;
		double MAX_DISPLACEMENT = 0.0;
		for (int k=1; k<Nz-1; k++){
			for (int j=1; j<Ny-1; j++){
				for (int i=1; i<Nx-1; i++){
					double walldist=BoundaryDist(i,j,k);
					double wallweight = 1.0 / (1+exp(-5.f*(walldist-1.f))); 
					//wallweight = 1.0;
					if (fabs(wallweight*morph_delta) > MAX_DISPLACEMENT) MAX_DISPLACEMENT= fabs(wallweight*morph_delta);
					
					if (Dist(i,j,k) - wallweight*morph_delta < 0.0){
						count+=1.0;
					}
				}
			}
		}
		count = Dm->Comm.sumReduce( count );
		MAX_DISPLACEMENT = Dm->Comm.maxReduce( MAX_DISPLACEMENT );
		GrowthEstimate = count - count_original;
		ERROR = fabs((GrowthEstimate-TargetGrowth) /TargetGrowth);

		if (rank == 0) printf("     delta=%f, growth=%f, max. displacement = %f \n",morph_delta, GrowthEstimate, MAX_DISPLACEMENT);
		// Now adjust morph_delta
		double step_size = (TargetGrowth - GrowthEstimate)*(morph_delta - morph_delta_previous) / (GrowthEstimate - GrowthPrevious);
		GrowthPrevious = GrowthEstimate;
		morph_delta_previous = morph_delta;
		morph_delta += step_size;
		if (morph_delta / morph_delta_previous > 2.0 ) morph_delta = morph_delta_previous*2.0;
				
		//MAX_DISPLACEMENT *= max(TargetGrowth/GrowthEstimate,1.25);
		
		if (morph_delta > 0.0 ){
			// object is growing
			if (MAX_DISPLACEMENT > 3.0 ){
				morph_delta = 3.0;
				COUNT_FOR_LOOP = 100; // exit loop if displacement is too large
			}
		}
		else{
			// object is shrinking
			if (MAX_DISPLACEMENT > 1.0 ){
				morph_delta = -1.0;
				COUNT_FOR_LOOP = 100; // exit loop if displacement is too large
			}
		}
	}
	if (rank == 0) printf("Final delta=%f \n",morph_delta);

	count = 0.0;
	for (int k=1; k<Nz-1; k++){
		for (int j=1; j<Ny-1; j++){
			for (int i=1; i<Nx-1; i++){
				double walldist=BoundaryDist(i,j,k);
				double wallweight = 1.0 / (1+exp(-5.f*(walldist-1.f))); 
				//wallweight = 1.0;
				Dist(i,j,k) -= wallweight*morph_delta;
				if (Dist(i,j,k) < 0.0)	count+=1.0;
			}
		}
	}
	count = Dm->Comm.sumReduce( count );

	return count;
}


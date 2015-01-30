// Header file for two-phase averaging class
#include "pmmc.h"
#include "Domain.h"
#include "Communication.h"

class TwoPhase{

	//...........................................................................
	int n_nw_pts,n_ns_pts,n_ws_pts,n_nws_pts,n_local_sol_pts,n_local_nws_pts;
	int n_nw_tris,n_ns_tris,n_ws_tris,n_nws_seg,n_local_sol_tris;
	int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};
	//...........................................................................
	int nc;
	int kstart,kfinish;

	double fluid_isovalue, solid_isovalue;
	double Volume;
	// initialize lists for vertices for surfaces, common line
	DTMutableList<Point> nw_pts;
	DTMutableList<Point> ns_pts;
	DTMutableList<Point> ws_pts;
	DTMutableList<Point> nws_pts;
	DTMutableList<Point> local_sol_pts;
	DTMutableList<Point> local_nws_pts;
	DTMutableList<Point> tmp;

	// initialize triangle lists for surfaces
	IntArray nw_tris;
	IntArray ns_tris;
	IntArray ws_tris;
	IntArray nws_seg;
	IntArray local_sol_tris;

	IntArray cubeList;

	// Temporary storage arrays
	DoubleArray CubeValues;
	DoubleArray Values;
	DoubleArray DistanceValues;
	DoubleArray KGwns_values;
	DoubleArray KNwns_values;
	DoubleArray InterfaceSpeed;
	DoubleArray NormalVector;

	// CSV / text file where time history of averages is saved
	FILE *TIMELOG;

	// Buffers for communication
	// Fill in the phase MeshData from neighboring processors
	// Note that the communicator is associated with the domain Dm
	double *sendMeshData_x, *sendMeshData_y, *sendMeshData_z, *sendMeshData_X, *sendMeshData_Y, *sendMeshData_Z;
	double *sendMeshData_xy, *sendMeshData_yz, *sendMeshData_xz, *sendMeshData_Xy, *sendMeshData_Yz, *sendMeshData_xZ;
	double *sendMeshData_xY, *sendMeshData_yZ, *sendMeshData_Xz, *sendMeshData_XY, *sendMeshData_YZ, *sendMeshData_XZ;
	double *recvMeshData_x, *recvMeshData_y, *recvMeshData_z, *recvMeshData_X, *recvMeshData_Y, *recvMeshData_Z;
	double *recvMeshData_xy, *recvMeshData_yz, *recvMeshData_xz, *recvMeshData_Xy, *recvMeshData_Yz, *recvMeshData_xZ;
	double *recvMeshData_xY, *recvMeshData_yZ, *recvMeshData_Xz, *recvMeshData_XY, *recvMeshData_YZ, *recvMeshData_XZ;

public:
	Domain& Dm;
	int ncubes;
	//...........................................................................
	// Averaging variables
	//...........................................................................
	// local averages (to each MPI process)
	double trimdist; 						// pixel distance to trim surface for specified averages
	double porosity,poreVol;
	double awn,ans,aws,lwns;
	double wp_volume,nwp_volume;
	double As, dummy;
	double vol_w, vol_n;						// volumes the exclude the interfacial region
	double sat_w, sat_w_previous;
	double pan,paw;								// local phase averaged pressure
	// Global averages (all processes)
	double pan_global,paw_global;				// local phase averaged pressure
	double vol_w_global, vol_n_global;			// volumes the exclude the interfacial region
	double awn_global,ans_global,aws_global;
	double lwns_global;
	double efawns,efawns_global;				// averaged contact angle
	double Jwn,Jwn_global;						// average mean curavture - wn interface
	double Kwn,Kwn_global;						// average Gaussian curavture - wn interface
	double KNwns,KNwns_global;					// wns common curve normal curavture
	double KGwns,KGwns_global;					// wns common curve geodesic curavture
	double trawn,trawn_global;					// trimmed interfacial area
	double trJwn,trJwn_global;					// trimmed interfacial area
	double trRwn,trRwn_global;					// trimmed interfacial area
	double nwp_volume_global;					// volume for the non-wetting phase
	double wp_volume_global;					// volume for the wetting phase
	double As_global;
	double dEs,dAwn,dAns;						// Global surface energy (calculated by rank=0)
	DoubleArray van;
	DoubleArray vaw;
	DoubleArray vawn;
	DoubleArray vawns;
	DoubleArray Gwn;
	DoubleArray Gns;
	DoubleArray Gws;
	DoubleArray van_global;
	DoubleArray vaw_global;
	DoubleArray vawn_global;
	DoubleArray vawns_global;
	DoubleArray Gwn_global;
	DoubleArray Gns_global;
	DoubleArray Gws_global;
	//...........................................................................
	//...........................................................................
	int Nx,Ny,Nz;
	IntArray LocalBlobID;
	DoubleArray SDn;
	DoubleArray SDs;
	DoubleArray Phase;
	DoubleArray Press;
	DoubleArray dPdt;
	DoubleArray MeanCurvature;
	DoubleArray GaussCurvature;
	DoubleArray SDs_x;		// Gradient of the signed distance
	DoubleArray SDs_y;
	DoubleArray SDs_z;
	DoubleArray SDn_x;		// Gradient of the signed distance
	DoubleArray SDn_y;
	DoubleArray SDn_z;
	DoubleArray DelPhi;		// Magnitude of Gradient of the phase indicator field
	DoubleArray Phase_tplus;
	DoubleArray Phase_tminus;
	DoubleArray Vel_x;		// Velocity
	DoubleArray Vel_y;
	DoubleArray Vel_z;
	//...........................................................................
	TwoPhase(Domain &dm) : Dm(dm){
		Nx=dm.Nx; Ny=dm.Ny; Nz=dm.Nz;
		Volume=(Nx-2)*(Ny-2)*(Nz-2)*Dm.nprocx*Dm.nprocy*Dm.nprocz*1.0;

		ncubes=(Nx-2)*(Ny-2)*(Nz-2);
		cubeList.New(3,ncubes);

		// Global arrays
		LocalBlobID.New(Nx,Ny,Nz);
		SDn.New(Nx,Ny,Nz);
		SDs.New(Nx,Ny,Nz);
		Phase.New(Nx,Ny,Nz);
		Press.New(Nx,Ny,Nz);
		dPdt.New(Nx,Ny,Nz);
		MeanCurvature.New(Nx,Ny,Nz);
		GaussCurvature.New(Nx,Ny,Nz);
		SDs_x.New(Nx,Ny,Nz);		// Gradient of the signed distance
		SDs_y.New(Nx,Ny,Nz);
		SDs_z.New(Nx,Ny,Nz);
		SDn_x.New(Nx,Ny,Nz);		// Gradient of the signed distance
		SDn_y.New(Nx,Ny,Nz);
		SDn_z.New(Nx,Ny,Nz);
		DelPhi.New(Nx,Ny,Nz);
		Phase_tplus.New(Nx,Ny,Nz);
		Phase_tminus.New(Nx,Ny,Nz);
		Vel_x.New(Nx,Ny,Nz);		// Gradient of the phase indicator field
		Vel_y.New(Nx,Ny,Nz);
		Vel_z.New(Nx,Ny,Nz);
		//.........................................
		// Allocate cube storage space
		CubeValues.New(2,2,2);
		nw_tris.New(3,20);
		ns_tris.New(3,20);
		ws_tris.New(3,20);
		nws_seg.New(2,20);
		local_sol_tris.New(3,18);
		nw_pts=DTMutableList<Point>(20);
		ns_pts=DTMutableList<Point>(20);
		ws_pts=DTMutableList<Point>(20);
		nws_pts=DTMutableList<Point>(20);
		local_nws_pts=DTMutableList<Point>(20);
		local_sol_pts=DTMutableList<Point>(20);
		tmp=DTMutableList<Point>(20);
		//.........................................
		Values.New(20);
		DistanceValues.New(20);
		KGwns_values.New(20);
		KNwns_values.New(20);
		InterfaceSpeed.New(20);
		NormalVector.New(60);
		//.........................................
		van.New(3);
		vaw.New(3);
		vawn.New(3);
		vawns.New(3);
		Gwn.New(6);
		Gns.New(6);
		Gws.New(6);
		van_global.New(3);
		vaw_global.New(3);
		vawn_global.New(3);
		vawns_global.New(3);
		Gwn_global.New(6);
		Gns_global.New(6);
		Gws_global.New(6);
		//.........................................
		if (Dm.rank==0){
			TIMELOG = fopen("timelog.tcat","a+");
			if (fseek(TIMELOG,0,SEEK_SET) == fseek(TIMELOG,0,SEEK_CUR)){
				// If timelog is empty, write a short header to list the averages
				//fprintf(TIMELOG,"--------------------------------------------------------------------------------------\n");
				fprintf(TIMELOG,"time dEs ");								// Timestep, Change in Surface Energy
				fprintf(TIMELOG,"sw pw pn awn ans aws Jwn Kwn lwns sgkvpmawns KNwns KGwns ");	// Scalar averages
				fprintf(TIMELOG,"vawx vawy vawz vanx vany vanz ");			// Velocity averages
				fprintf(TIMELOG,"vawnx vawny vawnz vawnsx vawnsy vawnsz ");
				fprintf(TIMELOG,"Gwnxx Gwnyy Gwnzz Gwnxy Gwnxz Gwnyz ");				// Orientation tensors
				fprintf(TIMELOG,"Gwsxx Gwsyy Gwszz Gwsxy Gwsxz Gwsyz ");
				fprintf(TIMELOG,"Gnsxx Gnsyy Gnszz Gnsxy Gnsxz Gnsyz ");
				fprintf(TIMELOG,"trawn trJwn trRwn\n");								// trimmed curvature for wn surface
				//fprintf(TIMELOG,"--------------------------------------------------------------------------------------\n");
			}
		}
		//.........................................
		// Allocation buffers for mesh communication
		// send buffers
		sendMeshData_x = new double [Dm.sendCount_x];
		sendMeshData_y = new double [Dm.sendCount_y];
		sendMeshData_z = new double [Dm.sendCount_z];
		sendMeshData_X = new double [Dm.sendCount_X];
		sendMeshData_Y = new double [Dm.sendCount_Y];
		sendMeshData_Z = new double [Dm.sendCount_Z];
		sendMeshData_xy = new double [Dm.sendCount_xy];
		sendMeshData_yz = new double [Dm.sendCount_yz];
		sendMeshData_xz = new double [Dm.sendCount_xz];
		sendMeshData_Xy = new double [Dm.sendCount_Xy];
		sendMeshData_Yz = new double [Dm.sendCount_Yz];
		sendMeshData_xZ = new double [Dm.sendCount_xZ];
		sendMeshData_xY = new double [Dm.sendCount_xY];
		sendMeshData_yZ = new double [Dm.sendCount_yZ];
		sendMeshData_Xz = new double [Dm.sendCount_Xz];
		sendMeshData_XY = new double [Dm.sendCount_XY];
		sendMeshData_YZ = new double [Dm.sendCount_YZ];
		sendMeshData_XZ = new double [Dm.sendCount_XZ];
		//......................................................................................
		// recv buffers
		recvMeshData_x = new double [Dm.recvCount_x];
		recvMeshData_y = new double [Dm.recvCount_y];
		recvMeshData_z = new double [Dm.recvCount_z];
		recvMeshData_X = new double [Dm.recvCount_X];
		recvMeshData_Y = new double [Dm.recvCount_Y];
		recvMeshData_Z = new double [Dm.recvCount_Z];
		recvMeshData_xy = new double [Dm.recvCount_xy];
		recvMeshData_yz = new double [Dm.recvCount_yz];
		recvMeshData_xz = new double [Dm.recvCount_xz];
		recvMeshData_Xy = new double [Dm.recvCount_Xy];
		recvMeshData_xZ = new double [Dm.recvCount_xZ];
		recvMeshData_xY = new double [Dm.recvCount_xY];
		recvMeshData_yZ = new double [Dm.recvCount_yZ];
		recvMeshData_Yz = new double [Dm.recvCount_Yz];
		recvMeshData_Xz = new double [Dm.recvCount_Xz];
		recvMeshData_XY = new double [Dm.recvCount_XY];
		recvMeshData_YZ = new double [Dm.recvCount_YZ];
		recvMeshData_XZ = new double [Dm.recvCount_XZ];
		//......................................................................................

	}
	~TwoPhase(){
		delete sendMeshData_x;
		delete sendMeshData_y;
		delete sendMeshData_z;
		delete sendMeshData_X;
		delete sendMeshData_Y;
		delete sendMeshData_Z;
		delete sendMeshData_xy;
		delete sendMeshData_xY;
		delete sendMeshData_Xy;
		delete sendMeshData_XY;
		delete sendMeshData_xz;
		delete sendMeshData_xZ;
		delete sendMeshData_Xz;
		delete sendMeshData_XZ;
		delete sendMeshData_yz;
		delete sendMeshData_yZ;
		delete sendMeshData_Yz;
		delete sendMeshData_YZ;
		delete recvMeshData_x;
		delete recvMeshData_y;
		delete recvMeshData_z;
		delete recvMeshData_X;
		delete recvMeshData_Y;
		delete recvMeshData_Z;
		delete recvMeshData_xy;
		delete recvMeshData_xY;
		delete recvMeshData_Xy;
		delete recvMeshData_XY;
		delete recvMeshData_xz;
		delete recvMeshData_xZ;
		delete recvMeshData_Xz;
		delete recvMeshData_XZ;
		delete recvMeshData_yz;
		delete recvMeshData_yZ;
		delete recvMeshData_Yz;
		delete recvMeshData_YZ;
	}

	void Initialize();
	void SetupCubes(Domain &Dm);
	void UpdateMeshValues();
	void UpdateSolid();
	void ComputeDelPhi();
	void ColorToSignedDistance(double Beta, double *ColorData, double *DistData);
	void ComputeLocal();
	void Reduce();
	void NonDimensionalize(double D, double viscosity, double IFT);
	void PrintAll(int timestep);

};

void TwoPhase::ColorToSignedDistance(double Beta, double *ColorData, double *DistData){

	double temp=0.5/Beta;
/*	for (int n=0; n<Nx*Ny*Nz; n++){
		double value = ColorData[n];
		DistData[n] = temp*log((1.0+value)/(1.0-value));
	}
*/
	for (int n=0; n<Nx*Ny*Nz; n++)	DistData[n] = ColorData[n];
}

void TwoPhase::ComputeDelPhi(){

	int i,j,k;
	double fx,fy,fz;

	CommunicateMeshHalo(Phase, Dm.Comm,
			sendMeshData_x,sendMeshData_y,sendMeshData_z,sendMeshData_X,sendMeshData_Y,sendMeshData_Z,
			sendMeshData_xy,sendMeshData_XY,sendMeshData_xY,sendMeshData_Xy,sendMeshData_xz,sendMeshData_XZ,
			sendMeshData_xZ,sendMeshData_Xz,sendMeshData_yz,sendMeshData_YZ,sendMeshData_yZ,sendMeshData_Yz,
			recvMeshData_x,recvMeshData_y,recvMeshData_z,recvMeshData_X,recvMeshData_Y,recvMeshData_Z,
			recvMeshData_xy,recvMeshData_XY,recvMeshData_xY,recvMeshData_Xy,recvMeshData_xz,recvMeshData_XZ,
			recvMeshData_xZ,recvMeshData_Xz,recvMeshData_yz,recvMeshData_YZ,recvMeshData_yZ,recvMeshData_Yz,
			Dm.sendList_x,Dm.sendList_y,Dm.sendList_z,Dm.sendList_X,Dm.sendList_Y,Dm.sendList_Z,
			Dm.sendList_xy,Dm.sendList_XY,Dm.sendList_xY,Dm.sendList_Xy,Dm.sendList_xz,Dm.sendList_XZ,
			Dm.sendList_xZ,Dm.sendList_Xz,Dm.sendList_yz,Dm.sendList_YZ,Dm.sendList_yZ,Dm.sendList_Yz,
			Dm.sendCount_x,Dm.sendCount_y,Dm.sendCount_z,Dm.sendCount_X,Dm.sendCount_Y,Dm.sendCount_Z,
			Dm.sendCount_xy,Dm.sendCount_XY,Dm.sendCount_xY,Dm.sendCount_Xy,Dm.sendCount_xz,Dm.sendCount_XZ,
			Dm.sendCount_xZ,Dm.sendCount_Xz,Dm.sendCount_yz,Dm.sendCount_YZ,Dm.sendCount_yZ,Dm.sendCount_Yz,
			Dm.recvList_x,Dm.recvList_y,Dm.recvList_z,Dm.recvList_X,Dm.recvList_Y,Dm.recvList_Z,
			Dm.recvList_xy,Dm.recvList_XY,Dm.recvList_xY,Dm.recvList_Xy,Dm.recvList_xz,Dm.recvList_XZ,
			Dm.recvList_xZ,Dm.recvList_Xz,Dm.recvList_yz,Dm.recvList_YZ,Dm.recvList_yZ,Dm.recvList_Yz,
			Dm.recvCount_x,Dm.recvCount_y,Dm.recvCount_z,Dm.recvCount_X,Dm.recvCount_Y,Dm.recvCount_Z,
			Dm.recvCount_xy,Dm.recvCount_XY,Dm.recvCount_xY,Dm.recvCount_Xy,Dm.recvCount_xz,Dm.recvCount_XZ,
			Dm.recvCount_xZ,Dm.recvCount_Xz,Dm.recvCount_yz,Dm.recvCount_YZ,Dm.recvCount_yZ,Dm.recvCount_Yz,
			Dm.rank_x,Dm.rank_y,Dm.rank_z,Dm.rank_X,Dm.rank_Y,Dm.rank_Z,Dm.rank_xy,Dm.rank_XY,Dm.rank_xY,
			Dm.rank_Xy,Dm.rank_xz,Dm.rank_XZ,Dm.rank_xZ,Dm.rank_Xz,Dm.rank_yz,Dm.rank_YZ,Dm.rank_yZ,Dm.rank_Yz);

	for (k=1; k<Nz-1; k++){
		for (j=1; j<Ny-1; j++){
			for (i=1; i<Nx-1; i++){
				// Compute all of the derivatives using finite differences
				fx = 0.5*(Phase(i+1,j,k) - Phase(i-1,j,k));
				fy = 0.5*(Phase(i,j+1,k) - Phase(i,j-1,k));
				fz = 0.5*(Phase(i,j,k+1) - Phase(i,j,k-1));
				DelPhi(i,j,k) = sqrt(fx*fx+fy*fy+fz*fz);
			}
		}
	}
}

void TwoPhase::Initialize(){

	trimdist=1.0;
	fluid_isovalue=solid_isovalue=0.0;
	// Initialize the averaged quantities
	awn = aws = ans = lwns = 0.0;
	nwp_volume = wp_volume = 0.0;
	As = 0.0;
	pan = paw = 0.0;
	vaw(0) = vaw(1) = vaw(2) = 0.0;
	van(0) = van(1) = van(2) = 0.0;
	vawn(0) = vawn(1) = vawn(2) = 0.0;
	vawns(0) = vawns(1) = vawns(2) = 0.0;
	Gwn(0) = Gwn(1) = Gwn(2) = 0.0;
	Gwn(3) = Gwn(4) = Gwn(5) = 0.0;
	Gws(0) = Gws(1) = Gws(2) = 0.0;
	Gws(3) = Gws(4) = Gws(5) = 0.0;
	Gns(0) = Gns(1) = Gns(2) = 0.0;
	Gns(3) = Gns(4) = Gns(5) = 0.0;
	vol_w = vol_n =0.0;
	KGwns = KNwns = 0.0;
	Jwn = Kwn = efawns = 0.0;
	trJwn = trawn = trRwn = 0.0;
}

void TwoPhase::SetupCubes(Domain &Dm){
	int i,j,k;
	kstart = 1;
	kfinish = Nz-1;
	if (Dm.BoundaryCondition !=0 && Dm.kproc==0)			kstart = 4;
	if (Dm.BoundaryCondition !=0 && Dm.kproc==Dm.nprocz-1)	kfinish = Nz-4;
	nc=0;
	for (k=kstart; k<kfinish; k++){
		for (j=1; j<Ny-1; j++){
			for (i=1; i<Nx-1; i++){
				cubeList(0,nc) = i;
				cubeList(1,nc) = j;
				cubeList(2,nc) = k;
				nc++;
			}
		}
	}
	ncubes = nc;
}

void TwoPhase::UpdateSolid(){
	//...........................................................................
	// Gradient of the Signed Distance function
	//...........................................................................
	pmmc_MeshGradient(SDs,SDs_x,SDs_y,SDs_z,Nx,Ny,Nz);
	//...........................................................................
	CommunicateMeshHalo(SDs_x, Dm.Comm,
			sendMeshData_x,sendMeshData_y,sendMeshData_z,sendMeshData_X,sendMeshData_Y,sendMeshData_Z,
			sendMeshData_xy,sendMeshData_XY,sendMeshData_xY,sendMeshData_Xy,sendMeshData_xz,sendMeshData_XZ,
			sendMeshData_xZ,sendMeshData_Xz,sendMeshData_yz,sendMeshData_YZ,sendMeshData_yZ,sendMeshData_Yz,
			recvMeshData_x,recvMeshData_y,recvMeshData_z,recvMeshData_X,recvMeshData_Y,recvMeshData_Z,
			recvMeshData_xy,recvMeshData_XY,recvMeshData_xY,recvMeshData_Xy,recvMeshData_xz,recvMeshData_XZ,
			recvMeshData_xZ,recvMeshData_Xz,recvMeshData_yz,recvMeshData_YZ,recvMeshData_yZ,recvMeshData_Yz,
			Dm.sendList_x,Dm.sendList_y,Dm.sendList_z,Dm.sendList_X,Dm.sendList_Y,Dm.sendList_Z,
			Dm.sendList_xy,Dm.sendList_XY,Dm.sendList_xY,Dm.sendList_Xy,Dm.sendList_xz,Dm.sendList_XZ,
			Dm.sendList_xZ,Dm.sendList_Xz,Dm.sendList_yz,Dm.sendList_YZ,Dm.sendList_yZ,Dm.sendList_Yz,
			Dm.sendCount_x,Dm.sendCount_y,Dm.sendCount_z,Dm.sendCount_X,Dm.sendCount_Y,Dm.sendCount_Z,
			Dm.sendCount_xy,Dm.sendCount_XY,Dm.sendCount_xY,Dm.sendCount_Xy,Dm.sendCount_xz,Dm.sendCount_XZ,
			Dm.sendCount_xZ,Dm.sendCount_Xz,Dm.sendCount_yz,Dm.sendCount_YZ,Dm.sendCount_yZ,Dm.sendCount_Yz,
			Dm.recvList_x,Dm.recvList_y,Dm.recvList_z,Dm.recvList_X,Dm.recvList_Y,Dm.recvList_Z,
			Dm.recvList_xy,Dm.recvList_XY,Dm.recvList_xY,Dm.recvList_Xy,Dm.recvList_xz,Dm.recvList_XZ,
			Dm.recvList_xZ,Dm.recvList_Xz,Dm.recvList_yz,Dm.recvList_YZ,Dm.recvList_yZ,Dm.recvList_Yz,
			Dm.recvCount_x,Dm.recvCount_y,Dm.recvCount_z,Dm.recvCount_X,Dm.recvCount_Y,Dm.recvCount_Z,
			Dm.recvCount_xy,Dm.recvCount_XY,Dm.recvCount_xY,Dm.recvCount_Xy,Dm.recvCount_xz,Dm.recvCount_XZ,
			Dm.recvCount_xZ,Dm.recvCount_Xz,Dm.recvCount_yz,Dm.recvCount_YZ,Dm.recvCount_yZ,Dm.recvCount_Yz,
			Dm.rank_x,Dm.rank_y,Dm.rank_z,Dm.rank_X,Dm.rank_Y,Dm.rank_Z,Dm.rank_xy,Dm.rank_XY,Dm.rank_xY,
			Dm.rank_Xy,Dm.rank_xz,Dm.rank_XZ,Dm.rank_xZ,Dm.rank_Xz,Dm.rank_yz,Dm.rank_YZ,Dm.rank_yZ,Dm.rank_Yz);
	//...........................................................................
	CommunicateMeshHalo(SDs_y, Dm.Comm,
			sendMeshData_x,sendMeshData_y,sendMeshData_z,sendMeshData_X,sendMeshData_Y,sendMeshData_Z,
			sendMeshData_xy,sendMeshData_XY,sendMeshData_xY,sendMeshData_Xy,sendMeshData_xz,sendMeshData_XZ,
			sendMeshData_xZ,sendMeshData_Xz,sendMeshData_yz,sendMeshData_YZ,sendMeshData_yZ,sendMeshData_Yz,
			recvMeshData_x,recvMeshData_y,recvMeshData_z,recvMeshData_X,recvMeshData_Y,recvMeshData_Z,
			recvMeshData_xy,recvMeshData_XY,recvMeshData_xY,recvMeshData_Xy,recvMeshData_xz,recvMeshData_XZ,
			recvMeshData_xZ,recvMeshData_Xz,recvMeshData_yz,recvMeshData_YZ,recvMeshData_yZ,recvMeshData_Yz,
			Dm.sendList_x,Dm.sendList_y,Dm.sendList_z,Dm.sendList_X,Dm.sendList_Y,Dm.sendList_Z,
			Dm.sendList_xy,Dm.sendList_XY,Dm.sendList_xY,Dm.sendList_Xy,Dm.sendList_xz,Dm.sendList_XZ,
			Dm.sendList_xZ,Dm.sendList_Xz,Dm.sendList_yz,Dm.sendList_YZ,Dm.sendList_yZ,Dm.sendList_Yz,
			Dm.sendCount_x,Dm.sendCount_y,Dm.sendCount_z,Dm.sendCount_X,Dm.sendCount_Y,Dm.sendCount_Z,
			Dm.sendCount_xy,Dm.sendCount_XY,Dm.sendCount_xY,Dm.sendCount_Xy,Dm.sendCount_xz,Dm.sendCount_XZ,
			Dm.sendCount_xZ,Dm.sendCount_Xz,Dm.sendCount_yz,Dm.sendCount_YZ,Dm.sendCount_yZ,Dm.sendCount_Yz,
			Dm.recvList_x,Dm.recvList_y,Dm.recvList_z,Dm.recvList_X,Dm.recvList_Y,Dm.recvList_Z,
			Dm.recvList_xy,Dm.recvList_XY,Dm.recvList_xY,Dm.recvList_Xy,Dm.recvList_xz,Dm.recvList_XZ,
			Dm.recvList_xZ,Dm.recvList_Xz,Dm.recvList_yz,Dm.recvList_YZ,Dm.recvList_yZ,Dm.recvList_Yz,
			Dm.recvCount_x,Dm.recvCount_y,Dm.recvCount_z,Dm.recvCount_X,Dm.recvCount_Y,Dm.recvCount_Z,
			Dm.recvCount_xy,Dm.recvCount_XY,Dm.recvCount_xY,Dm.recvCount_Xy,Dm.recvCount_xz,Dm.recvCount_XZ,
			Dm.recvCount_xZ,Dm.recvCount_Xz,Dm.recvCount_yz,Dm.recvCount_YZ,Dm.recvCount_yZ,Dm.recvCount_Yz,
			Dm.rank_x,Dm.rank_y,Dm.rank_z,Dm.rank_X,Dm.rank_Y,Dm.rank_Z,Dm.rank_xy,Dm.rank_XY,Dm.rank_xY,
			Dm.rank_Xy,Dm.rank_xz,Dm.rank_XZ,Dm.rank_xZ,Dm.rank_Xz,Dm.rank_yz,Dm.rank_YZ,Dm.rank_yZ,Dm.rank_Yz);
	//...........................................................................
	CommunicateMeshHalo(SDs_z, Dm.Comm,
			sendMeshData_x,sendMeshData_y,sendMeshData_z,sendMeshData_X,sendMeshData_Y,sendMeshData_Z,
			sendMeshData_xy,sendMeshData_XY,sendMeshData_xY,sendMeshData_Xy,sendMeshData_xz,sendMeshData_XZ,
			sendMeshData_xZ,sendMeshData_Xz,sendMeshData_yz,sendMeshData_YZ,sendMeshData_yZ,sendMeshData_Yz,
			recvMeshData_x,recvMeshData_y,recvMeshData_z,recvMeshData_X,recvMeshData_Y,recvMeshData_Z,
			recvMeshData_xy,recvMeshData_XY,recvMeshData_xY,recvMeshData_Xy,recvMeshData_xz,recvMeshData_XZ,
			recvMeshData_xZ,recvMeshData_Xz,recvMeshData_yz,recvMeshData_YZ,recvMeshData_yZ,recvMeshData_Yz,
			Dm.sendList_x,Dm.sendList_y,Dm.sendList_z,Dm.sendList_X,Dm.sendList_Y,Dm.sendList_Z,
			Dm.sendList_xy,Dm.sendList_XY,Dm.sendList_xY,Dm.sendList_Xy,Dm.sendList_xz,Dm.sendList_XZ,
			Dm.sendList_xZ,Dm.sendList_Xz,Dm.sendList_yz,Dm.sendList_YZ,Dm.sendList_yZ,Dm.sendList_Yz,
			Dm.sendCount_x,Dm.sendCount_y,Dm.sendCount_z,Dm.sendCount_X,Dm.sendCount_Y,Dm.sendCount_Z,
			Dm.sendCount_xy,Dm.sendCount_XY,Dm.sendCount_xY,Dm.sendCount_Xy,Dm.sendCount_xz,Dm.sendCount_XZ,
			Dm.sendCount_xZ,Dm.sendCount_Xz,Dm.sendCount_yz,Dm.sendCount_YZ,Dm.sendCount_yZ,Dm.sendCount_Yz,
			Dm.recvList_x,Dm.recvList_y,Dm.recvList_z,Dm.recvList_X,Dm.recvList_Y,Dm.recvList_Z,
			Dm.recvList_xy,Dm.recvList_XY,Dm.recvList_xY,Dm.recvList_Xy,Dm.recvList_xz,Dm.recvList_XZ,
			Dm.recvList_xZ,Dm.recvList_Xz,Dm.recvList_yz,Dm.recvList_YZ,Dm.recvList_yZ,Dm.recvList_Yz,
			Dm.recvCount_x,Dm.recvCount_y,Dm.recvCount_z,Dm.recvCount_X,Dm.recvCount_Y,Dm.recvCount_Z,
			Dm.recvCount_xy,Dm.recvCount_XY,Dm.recvCount_xY,Dm.recvCount_Xy,Dm.recvCount_xz,Dm.recvCount_XZ,
			Dm.recvCount_xZ,Dm.recvCount_Xz,Dm.recvCount_yz,Dm.recvCount_YZ,Dm.recvCount_yZ,Dm.recvCount_Yz,
			Dm.rank_x,Dm.rank_y,Dm.rank_z,Dm.rank_X,Dm.rank_Y,Dm.rank_Z,Dm.rank_xy,Dm.rank_XY,Dm.rank_xY,
			Dm.rank_Xy,Dm.rank_xz,Dm.rank_XZ,Dm.rank_xZ,Dm.rank_Xz,Dm.rank_yz,Dm.rank_YZ,Dm.rank_yZ,Dm.rank_Yz);
	//...........................................................................
}

void TwoPhase::UpdateMeshValues(){

	//...........................................................................
	// Compute the gradients of the phase indicator and signed distance fields
	pmmc_MeshGradient(SDs,SDs_x,SDs_y,SDs_z,Nx,Ny,Nz);
	pmmc_MeshGradient(SDn,SDn_x,SDn_y,SDn_z,Nx,Ny,Nz);
	//...........................................................................
	// Compute the mesh curvature of the phase indicator field
	pmmc_MeshCurvature(SDn, MeanCurvature, GaussCurvature, Nx, Ny, Nz);
	//...........................................................................
	// Update the time derivative of non-dimensional density field
	// Map Phase_tplus and Phase_tminus
	for (int n=0; n<Nx*Ny*Nz; n++)	dPdt(n) = 0.1*(Phase_tplus(n) - Phase_tminus(n));
	//...........................................................................
	CommunicateMeshHalo(Vel_x, Dm.Comm,
			sendMeshData_x,sendMeshData_y,sendMeshData_z,sendMeshData_X,sendMeshData_Y,sendMeshData_Z,
			sendMeshData_xy,sendMeshData_XY,sendMeshData_xY,sendMeshData_Xy,sendMeshData_xz,sendMeshData_XZ,
			sendMeshData_xZ,sendMeshData_Xz,sendMeshData_yz,sendMeshData_YZ,sendMeshData_yZ,sendMeshData_Yz,
			recvMeshData_x,recvMeshData_y,recvMeshData_z,recvMeshData_X,recvMeshData_Y,recvMeshData_Z,
			recvMeshData_xy,recvMeshData_XY,recvMeshData_xY,recvMeshData_Xy,recvMeshData_xz,recvMeshData_XZ,
			recvMeshData_xZ,recvMeshData_Xz,recvMeshData_yz,recvMeshData_YZ,recvMeshData_yZ,recvMeshData_Yz,
			Dm.sendList_x,Dm.sendList_y,Dm.sendList_z,Dm.sendList_X,Dm.sendList_Y,Dm.sendList_Z,
			Dm.sendList_xy,Dm.sendList_XY,Dm.sendList_xY,Dm.sendList_Xy,Dm.sendList_xz,Dm.sendList_XZ,
			Dm.sendList_xZ,Dm.sendList_Xz,Dm.sendList_yz,Dm.sendList_YZ,Dm.sendList_yZ,Dm.sendList_Yz,
			Dm.sendCount_x,Dm.sendCount_y,Dm.sendCount_z,Dm.sendCount_X,Dm.sendCount_Y,Dm.sendCount_Z,
			Dm.sendCount_xy,Dm.sendCount_XY,Dm.sendCount_xY,Dm.sendCount_Xy,Dm.sendCount_xz,Dm.sendCount_XZ,
			Dm.sendCount_xZ,Dm.sendCount_Xz,Dm.sendCount_yz,Dm.sendCount_YZ,Dm.sendCount_yZ,Dm.sendCount_Yz,
			Dm.recvList_x,Dm.recvList_y,Dm.recvList_z,Dm.recvList_X,Dm.recvList_Y,Dm.recvList_Z,
			Dm.recvList_xy,Dm.recvList_XY,Dm.recvList_xY,Dm.recvList_Xy,Dm.recvList_xz,Dm.recvList_XZ,
			Dm.recvList_xZ,Dm.recvList_Xz,Dm.recvList_yz,Dm.recvList_YZ,Dm.recvList_yZ,Dm.recvList_Yz,
			Dm.recvCount_x,Dm.recvCount_y,Dm.recvCount_z,Dm.recvCount_X,Dm.recvCount_Y,Dm.recvCount_Z,
			Dm.recvCount_xy,Dm.recvCount_XY,Dm.recvCount_xY,Dm.recvCount_Xy,Dm.recvCount_xz,Dm.recvCount_XZ,
			Dm.recvCount_xZ,Dm.recvCount_Xz,Dm.recvCount_yz,Dm.recvCount_YZ,Dm.recvCount_yZ,Dm.recvCount_Yz,
			Dm.rank_x,Dm.rank_y,Dm.rank_z,Dm.rank_X,Dm.rank_Y,Dm.rank_Z,Dm.rank_xy,Dm.rank_XY,Dm.rank_xY,
			Dm.rank_Xy,Dm.rank_xz,Dm.rank_XZ,Dm.rank_xZ,Dm.rank_Xz,Dm.rank_yz,Dm.rank_YZ,Dm.rank_yZ,Dm.rank_Yz);
	//...........................................................................
	// Velocity
	//...........................................................................
	CommunicateMeshHalo(Press, Dm.Comm,
			sendMeshData_x,sendMeshData_y,sendMeshData_z,sendMeshData_X,sendMeshData_Y,sendMeshData_Z,
			sendMeshData_xy,sendMeshData_XY,sendMeshData_xY,sendMeshData_Xy,sendMeshData_xz,sendMeshData_XZ,
			sendMeshData_xZ,sendMeshData_Xz,sendMeshData_yz,sendMeshData_YZ,sendMeshData_yZ,sendMeshData_Yz,
			recvMeshData_x,recvMeshData_y,recvMeshData_z,recvMeshData_X,recvMeshData_Y,recvMeshData_Z,
			recvMeshData_xy,recvMeshData_XY,recvMeshData_xY,recvMeshData_Xy,recvMeshData_xz,recvMeshData_XZ,
			recvMeshData_xZ,recvMeshData_Xz,recvMeshData_yz,recvMeshData_YZ,recvMeshData_yZ,recvMeshData_Yz,
			Dm.sendList_x,Dm.sendList_y,Dm.sendList_z,Dm.sendList_X,Dm.sendList_Y,Dm.sendList_Z,
			Dm.sendList_xy,Dm.sendList_XY,Dm.sendList_xY,Dm.sendList_Xy,Dm.sendList_xz,Dm.sendList_XZ,
			Dm.sendList_xZ,Dm.sendList_Xz,Dm.sendList_yz,Dm.sendList_YZ,Dm.sendList_yZ,Dm.sendList_Yz,
			Dm.sendCount_x,Dm.sendCount_y,Dm.sendCount_z,Dm.sendCount_X,Dm.sendCount_Y,Dm.sendCount_Z,
			Dm.sendCount_xy,Dm.sendCount_XY,Dm.sendCount_xY,Dm.sendCount_Xy,Dm.sendCount_xz,Dm.sendCount_XZ,
			Dm.sendCount_xZ,Dm.sendCount_Xz,Dm.sendCount_yz,Dm.sendCount_YZ,Dm.sendCount_yZ,Dm.sendCount_Yz,
			Dm.recvList_x,Dm.recvList_y,Dm.recvList_z,Dm.recvList_X,Dm.recvList_Y,Dm.recvList_Z,
			Dm.recvList_xy,Dm.recvList_XY,Dm.recvList_xY,Dm.recvList_Xy,Dm.recvList_xz,Dm.recvList_XZ,
			Dm.recvList_xZ,Dm.recvList_Xz,Dm.recvList_yz,Dm.recvList_YZ,Dm.recvList_yZ,Dm.recvList_Yz,
			Dm.recvCount_x,Dm.recvCount_y,Dm.recvCount_z,Dm.recvCount_X,Dm.recvCount_Y,Dm.recvCount_Z,
			Dm.recvCount_xy,Dm.recvCount_XY,Dm.recvCount_xY,Dm.recvCount_Xy,Dm.recvCount_xz,Dm.recvCount_XZ,
			Dm.recvCount_xZ,Dm.recvCount_Xz,Dm.recvCount_yz,Dm.recvCount_YZ,Dm.recvCount_yZ,Dm.recvCount_Yz,
			Dm.rank_x,Dm.rank_y,Dm.rank_z,Dm.rank_X,Dm.rank_Y,Dm.rank_Z,Dm.rank_xy,Dm.rank_XY,Dm.rank_xY,
			Dm.rank_Xy,Dm.rank_xz,Dm.rank_XZ,Dm.rank_xZ,Dm.rank_Xz,Dm.rank_yz,Dm.rank_YZ,Dm.rank_yZ,Dm.rank_Yz);
	//...........................................................................
	//...........................................................................
	CommunicateMeshHalo(Vel_y, Dm.Comm,
			sendMeshData_x,sendMeshData_y,sendMeshData_z,sendMeshData_X,sendMeshData_Y,sendMeshData_Z,
			sendMeshData_xy,sendMeshData_XY,sendMeshData_xY,sendMeshData_Xy,sendMeshData_xz,sendMeshData_XZ,
			sendMeshData_xZ,sendMeshData_Xz,sendMeshData_yz,sendMeshData_YZ,sendMeshData_yZ,sendMeshData_Yz,
			recvMeshData_x,recvMeshData_y,recvMeshData_z,recvMeshData_X,recvMeshData_Y,recvMeshData_Z,
			recvMeshData_xy,recvMeshData_XY,recvMeshData_xY,recvMeshData_Xy,recvMeshData_xz,recvMeshData_XZ,
			recvMeshData_xZ,recvMeshData_Xz,recvMeshData_yz,recvMeshData_YZ,recvMeshData_yZ,recvMeshData_Yz,
			Dm.sendList_x,Dm.sendList_y,Dm.sendList_z,Dm.sendList_X,Dm.sendList_Y,Dm.sendList_Z,
			Dm.sendList_xy,Dm.sendList_XY,Dm.sendList_xY,Dm.sendList_Xy,Dm.sendList_xz,Dm.sendList_XZ,
			Dm.sendList_xZ,Dm.sendList_Xz,Dm.sendList_yz,Dm.sendList_YZ,Dm.sendList_yZ,Dm.sendList_Yz,
			Dm.sendCount_x,Dm.sendCount_y,Dm.sendCount_z,Dm.sendCount_X,Dm.sendCount_Y,Dm.sendCount_Z,
			Dm.sendCount_xy,Dm.sendCount_XY,Dm.sendCount_xY,Dm.sendCount_Xy,Dm.sendCount_xz,Dm.sendCount_XZ,
			Dm.sendCount_xZ,Dm.sendCount_Xz,Dm.sendCount_yz,Dm.sendCount_YZ,Dm.sendCount_yZ,Dm.sendCount_Yz,
			Dm.recvList_x,Dm.recvList_y,Dm.recvList_z,Dm.recvList_X,Dm.recvList_Y,Dm.recvList_Z,
			Dm.recvList_xy,Dm.recvList_XY,Dm.recvList_xY,Dm.recvList_Xy,Dm.recvList_xz,Dm.recvList_XZ,
			Dm.recvList_xZ,Dm.recvList_Xz,Dm.recvList_yz,Dm.recvList_YZ,Dm.recvList_yZ,Dm.recvList_Yz,
			Dm.recvCount_x,Dm.recvCount_y,Dm.recvCount_z,Dm.recvCount_X,Dm.recvCount_Y,Dm.recvCount_Z,
			Dm.recvCount_xy,Dm.recvCount_XY,Dm.recvCount_xY,Dm.recvCount_Xy,Dm.recvCount_xz,Dm.recvCount_XZ,
			Dm.recvCount_xZ,Dm.recvCount_Xz,Dm.recvCount_yz,Dm.recvCount_YZ,Dm.recvCount_yZ,Dm.recvCount_Yz,
			Dm.rank_x,Dm.rank_y,Dm.rank_z,Dm.rank_X,Dm.rank_Y,Dm.rank_Z,Dm.rank_xy,Dm.rank_XY,Dm.rank_xY,
			Dm.rank_Xy,Dm.rank_xz,Dm.rank_XZ,Dm.rank_xZ,Dm.rank_Xz,Dm.rank_yz,Dm.rank_YZ,Dm.rank_yZ,Dm.rank_Yz);
	//...........................................................................
	//...........................................................................
	CommunicateMeshHalo(Vel_z, Dm.Comm,
			sendMeshData_x,sendMeshData_y,sendMeshData_z,sendMeshData_X,sendMeshData_Y,sendMeshData_Z,
			sendMeshData_xy,sendMeshData_XY,sendMeshData_xY,sendMeshData_Xy,sendMeshData_xz,sendMeshData_XZ,
			sendMeshData_xZ,sendMeshData_Xz,sendMeshData_yz,sendMeshData_YZ,sendMeshData_yZ,sendMeshData_Yz,
			recvMeshData_x,recvMeshData_y,recvMeshData_z,recvMeshData_X,recvMeshData_Y,recvMeshData_Z,
			recvMeshData_xy,recvMeshData_XY,recvMeshData_xY,recvMeshData_Xy,recvMeshData_xz,recvMeshData_XZ,
			recvMeshData_xZ,recvMeshData_Xz,recvMeshData_yz,recvMeshData_YZ,recvMeshData_yZ,recvMeshData_Yz,
			Dm.sendList_x,Dm.sendList_y,Dm.sendList_z,Dm.sendList_X,Dm.sendList_Y,Dm.sendList_Z,
			Dm.sendList_xy,Dm.sendList_XY,Dm.sendList_xY,Dm.sendList_Xy,Dm.sendList_xz,Dm.sendList_XZ,
			Dm.sendList_xZ,Dm.sendList_Xz,Dm.sendList_yz,Dm.sendList_YZ,Dm.sendList_yZ,Dm.sendList_Yz,
			Dm.sendCount_x,Dm.sendCount_y,Dm.sendCount_z,Dm.sendCount_X,Dm.sendCount_Y,Dm.sendCount_Z,
			Dm.sendCount_xy,Dm.sendCount_XY,Dm.sendCount_xY,Dm.sendCount_Xy,Dm.sendCount_xz,Dm.sendCount_XZ,
			Dm.sendCount_xZ,Dm.sendCount_Xz,Dm.sendCount_yz,Dm.sendCount_YZ,Dm.sendCount_yZ,Dm.sendCount_Yz,
			Dm.recvList_x,Dm.recvList_y,Dm.recvList_z,Dm.recvList_X,Dm.recvList_Y,Dm.recvList_Z,
			Dm.recvList_xy,Dm.recvList_XY,Dm.recvList_xY,Dm.recvList_Xy,Dm.recvList_xz,Dm.recvList_XZ,
			Dm.recvList_xZ,Dm.recvList_Xz,Dm.recvList_yz,Dm.recvList_YZ,Dm.recvList_yZ,Dm.recvList_Yz,
			Dm.recvCount_x,Dm.recvCount_y,Dm.recvCount_z,Dm.recvCount_X,Dm.recvCount_Y,Dm.recvCount_Z,
			Dm.recvCount_xy,Dm.recvCount_XY,Dm.recvCount_xY,Dm.recvCount_Xy,Dm.recvCount_xz,Dm.recvCount_XZ,
			Dm.recvCount_xZ,Dm.recvCount_Xz,Dm.recvCount_yz,Dm.recvCount_YZ,Dm.recvCount_yZ,Dm.recvCount_Yz,
			Dm.rank_x,Dm.rank_y,Dm.rank_z,Dm.rank_X,Dm.rank_Y,Dm.rank_Z,Dm.rank_xy,Dm.rank_XY,Dm.rank_xY,
			Dm.rank_Xy,Dm.rank_xz,Dm.rank_XZ,Dm.rank_xZ,Dm.rank_Xz,Dm.rank_yz,Dm.rank_YZ,Dm.rank_yZ,Dm.rank_Yz);
	//...........................................................................
	// Mean Curvature
	//...........................................................................
	CommunicateMeshHalo(MeanCurvature, Dm.Comm,
			sendMeshData_x,sendMeshData_y,sendMeshData_z,sendMeshData_X,sendMeshData_Y,sendMeshData_Z,
			sendMeshData_xy,sendMeshData_XY,sendMeshData_xY,sendMeshData_Xy,sendMeshData_xz,sendMeshData_XZ,
			sendMeshData_xZ,sendMeshData_Xz,sendMeshData_yz,sendMeshData_YZ,sendMeshData_yZ,sendMeshData_Yz,
			recvMeshData_x,recvMeshData_y,recvMeshData_z,recvMeshData_X,recvMeshData_Y,recvMeshData_Z,
			recvMeshData_xy,recvMeshData_XY,recvMeshData_xY,recvMeshData_Xy,recvMeshData_xz,recvMeshData_XZ,
			recvMeshData_xZ,recvMeshData_Xz,recvMeshData_yz,recvMeshData_YZ,recvMeshData_yZ,recvMeshData_Yz,
			Dm.sendList_x,Dm.sendList_y,Dm.sendList_z,Dm.sendList_X,Dm.sendList_Y,Dm.sendList_Z,
			Dm.sendList_xy,Dm.sendList_XY,Dm.sendList_xY,Dm.sendList_Xy,Dm.sendList_xz,Dm.sendList_XZ,
			Dm.sendList_xZ,Dm.sendList_Xz,Dm.sendList_yz,Dm.sendList_YZ,Dm.sendList_yZ,Dm.sendList_Yz,
			Dm.sendCount_x,Dm.sendCount_y,Dm.sendCount_z,Dm.sendCount_X,Dm.sendCount_Y,Dm.sendCount_Z,
			Dm.sendCount_xy,Dm.sendCount_XY,Dm.sendCount_xY,Dm.sendCount_Xy,Dm.sendCount_xz,Dm.sendCount_XZ,
			Dm.sendCount_xZ,Dm.sendCount_Xz,Dm.sendCount_yz,Dm.sendCount_YZ,Dm.sendCount_yZ,Dm.sendCount_Yz,
			Dm.recvList_x,Dm.recvList_y,Dm.recvList_z,Dm.recvList_X,Dm.recvList_Y,Dm.recvList_Z,
			Dm.recvList_xy,Dm.recvList_XY,Dm.recvList_xY,Dm.recvList_Xy,Dm.recvList_xz,Dm.recvList_XZ,
			Dm.recvList_xZ,Dm.recvList_Xz,Dm.recvList_yz,Dm.recvList_YZ,Dm.recvList_yZ,Dm.recvList_Yz,
			Dm.recvCount_x,Dm.recvCount_y,Dm.recvCount_z,Dm.recvCount_X,Dm.recvCount_Y,Dm.recvCount_Z,
			Dm.recvCount_xy,Dm.recvCount_XY,Dm.recvCount_xY,Dm.recvCount_Xy,Dm.recvCount_xz,Dm.recvCount_XZ,
			Dm.recvCount_xZ,Dm.recvCount_Xz,Dm.recvCount_yz,Dm.recvCount_YZ,Dm.recvCount_yZ,Dm.recvCount_Yz,
			Dm.rank_x,Dm.rank_y,Dm.rank_z,Dm.rank_X,Dm.rank_Y,Dm.rank_Z,Dm.rank_xy,Dm.rank_XY,Dm.rank_xY,
			Dm.rank_Xy,Dm.rank_xz,Dm.rank_XZ,Dm.rank_xZ,Dm.rank_Xz,Dm.rank_yz,Dm.rank_YZ,Dm.rank_yZ,Dm.rank_Yz);
	//...........................................................................
	//...........................................................................
	// Gaussian Curvature
	//...........................................................................
	CommunicateMeshHalo(GaussCurvature, Dm.Comm,
			sendMeshData_x,sendMeshData_y,sendMeshData_z,sendMeshData_X,sendMeshData_Y,sendMeshData_Z,
			sendMeshData_xy,sendMeshData_XY,sendMeshData_xY,sendMeshData_Xy,sendMeshData_xz,sendMeshData_XZ,
			sendMeshData_xZ,sendMeshData_Xz,sendMeshData_yz,sendMeshData_YZ,sendMeshData_yZ,sendMeshData_Yz,
			recvMeshData_x,recvMeshData_y,recvMeshData_z,recvMeshData_X,recvMeshData_Y,recvMeshData_Z,
			recvMeshData_xy,recvMeshData_XY,recvMeshData_xY,recvMeshData_Xy,recvMeshData_xz,recvMeshData_XZ,
			recvMeshData_xZ,recvMeshData_Xz,recvMeshData_yz,recvMeshData_YZ,recvMeshData_yZ,recvMeshData_Yz,
			Dm.sendList_x,Dm.sendList_y,Dm.sendList_z,Dm.sendList_X,Dm.sendList_Y,Dm.sendList_Z,
			Dm.sendList_xy,Dm.sendList_XY,Dm.sendList_xY,Dm.sendList_Xy,Dm.sendList_xz,Dm.sendList_XZ,
			Dm.sendList_xZ,Dm.sendList_Xz,Dm.sendList_yz,Dm.sendList_YZ,Dm.sendList_yZ,Dm.sendList_Yz,
			Dm.sendCount_x,Dm.sendCount_y,Dm.sendCount_z,Dm.sendCount_X,Dm.sendCount_Y,Dm.sendCount_Z,
			Dm.sendCount_xy,Dm.sendCount_XY,Dm.sendCount_xY,Dm.sendCount_Xy,Dm.sendCount_xz,Dm.sendCount_XZ,
			Dm.sendCount_xZ,Dm.sendCount_Xz,Dm.sendCount_yz,Dm.sendCount_YZ,Dm.sendCount_yZ,Dm.sendCount_Yz,
			Dm.recvList_x,Dm.recvList_y,Dm.recvList_z,Dm.recvList_X,Dm.recvList_Y,Dm.recvList_Z,
			Dm.recvList_xy,Dm.recvList_XY,Dm.recvList_xY,Dm.recvList_Xy,Dm.recvList_xz,Dm.recvList_XZ,
			Dm.recvList_xZ,Dm.recvList_Xz,Dm.recvList_yz,Dm.recvList_YZ,Dm.recvList_yZ,Dm.recvList_Yz,
			Dm.recvCount_x,Dm.recvCount_y,Dm.recvCount_z,Dm.recvCount_X,Dm.recvCount_Y,Dm.recvCount_Z,
			Dm.recvCount_xy,Dm.recvCount_XY,Dm.recvCount_xY,Dm.recvCount_Xy,Dm.recvCount_xz,Dm.recvCount_XZ,
			Dm.recvCount_xZ,Dm.recvCount_Xz,Dm.recvCount_yz,Dm.recvCount_YZ,Dm.recvCount_yZ,Dm.recvCount_Yz,
			Dm.rank_x,Dm.rank_y,Dm.rank_z,Dm.rank_X,Dm.rank_Y,Dm.rank_Z,Dm.rank_xy,Dm.rank_XY,Dm.rank_xY,
			Dm.rank_Xy,Dm.rank_xz,Dm.rank_XZ,Dm.rank_xZ,Dm.rank_Xz,Dm.rank_yz,Dm.rank_YZ,Dm.rank_yZ,Dm.rank_Yz);
	//...........................................................................
	// Gradient of the phase indicator field
	//...........................................................................
	CommunicateMeshHalo(SDn_x, Dm.Comm,
			sendMeshData_x,sendMeshData_y,sendMeshData_z,sendMeshData_X,sendMeshData_Y,sendMeshData_Z,
			sendMeshData_xy,sendMeshData_XY,sendMeshData_xY,sendMeshData_Xy,sendMeshData_xz,sendMeshData_XZ,
			sendMeshData_xZ,sendMeshData_Xz,sendMeshData_yz,sendMeshData_YZ,sendMeshData_yZ,sendMeshData_Yz,
			recvMeshData_x,recvMeshData_y,recvMeshData_z,recvMeshData_X,recvMeshData_Y,recvMeshData_Z,
			recvMeshData_xy,recvMeshData_XY,recvMeshData_xY,recvMeshData_Xy,recvMeshData_xz,recvMeshData_XZ,
			recvMeshData_xZ,recvMeshData_Xz,recvMeshData_yz,recvMeshData_YZ,recvMeshData_yZ,recvMeshData_Yz,
			Dm.sendList_x,Dm.sendList_y,Dm.sendList_z,Dm.sendList_X,Dm.sendList_Y,Dm.sendList_Z,
			Dm.sendList_xy,Dm.sendList_XY,Dm.sendList_xY,Dm.sendList_Xy,Dm.sendList_xz,Dm.sendList_XZ,
			Dm.sendList_xZ,Dm.sendList_Xz,Dm.sendList_yz,Dm.sendList_YZ,Dm.sendList_yZ,Dm.sendList_Yz,
			Dm.sendCount_x,Dm.sendCount_y,Dm.sendCount_z,Dm.sendCount_X,Dm.sendCount_Y,Dm.sendCount_Z,
			Dm.sendCount_xy,Dm.sendCount_XY,Dm.sendCount_xY,Dm.sendCount_Xy,Dm.sendCount_xz,Dm.sendCount_XZ,
			Dm.sendCount_xZ,Dm.sendCount_Xz,Dm.sendCount_yz,Dm.sendCount_YZ,Dm.sendCount_yZ,Dm.sendCount_Yz,
			Dm.recvList_x,Dm.recvList_y,Dm.recvList_z,Dm.recvList_X,Dm.recvList_Y,Dm.recvList_Z,
			Dm.recvList_xy,Dm.recvList_XY,Dm.recvList_xY,Dm.recvList_Xy,Dm.recvList_xz,Dm.recvList_XZ,
			Dm.recvList_xZ,Dm.recvList_Xz,Dm.recvList_yz,Dm.recvList_YZ,Dm.recvList_yZ,Dm.recvList_Yz,
			Dm.recvCount_x,Dm.recvCount_y,Dm.recvCount_z,Dm.recvCount_X,Dm.recvCount_Y,Dm.recvCount_Z,
			Dm.recvCount_xy,Dm.recvCount_XY,Dm.recvCount_xY,Dm.recvCount_Xy,Dm.recvCount_xz,Dm.recvCount_XZ,
			Dm.recvCount_xZ,Dm.recvCount_Xz,Dm.recvCount_yz,Dm.recvCount_YZ,Dm.recvCount_yZ,Dm.recvCount_Yz,
			Dm.rank_x,Dm.rank_y,Dm.rank_z,Dm.rank_X,Dm.rank_Y,Dm.rank_Z,Dm.rank_xy,Dm.rank_XY,Dm.rank_xY,
			Dm.rank_Xy,Dm.rank_xz,Dm.rank_XZ,Dm.rank_xZ,Dm.rank_Xz,Dm.rank_yz,Dm.rank_YZ,Dm.rank_yZ,Dm.rank_Yz);
	//...........................................................................
	CommunicateMeshHalo(SDn_y, Dm.Comm,
			sendMeshData_x,sendMeshData_y,sendMeshData_z,sendMeshData_X,sendMeshData_Y,sendMeshData_Z,
			sendMeshData_xy,sendMeshData_XY,sendMeshData_xY,sendMeshData_Xy,sendMeshData_xz,sendMeshData_XZ,
			sendMeshData_xZ,sendMeshData_Xz,sendMeshData_yz,sendMeshData_YZ,sendMeshData_yZ,sendMeshData_Yz,
			recvMeshData_x,recvMeshData_y,recvMeshData_z,recvMeshData_X,recvMeshData_Y,recvMeshData_Z,
			recvMeshData_xy,recvMeshData_XY,recvMeshData_xY,recvMeshData_Xy,recvMeshData_xz,recvMeshData_XZ,
			recvMeshData_xZ,recvMeshData_Xz,recvMeshData_yz,recvMeshData_YZ,recvMeshData_yZ,recvMeshData_Yz,
			Dm.sendList_x,Dm.sendList_y,Dm.sendList_z,Dm.sendList_X,Dm.sendList_Y,Dm.sendList_Z,
			Dm.sendList_xy,Dm.sendList_XY,Dm.sendList_xY,Dm.sendList_Xy,Dm.sendList_xz,Dm.sendList_XZ,
			Dm.sendList_xZ,Dm.sendList_Xz,Dm.sendList_yz,Dm.sendList_YZ,Dm.sendList_yZ,Dm.sendList_Yz,
			Dm.sendCount_x,Dm.sendCount_y,Dm.sendCount_z,Dm.sendCount_X,Dm.sendCount_Y,Dm.sendCount_Z,
			Dm.sendCount_xy,Dm.sendCount_XY,Dm.sendCount_xY,Dm.sendCount_Xy,Dm.sendCount_xz,Dm.sendCount_XZ,
			Dm.sendCount_xZ,Dm.sendCount_Xz,Dm.sendCount_yz,Dm.sendCount_YZ,Dm.sendCount_yZ,Dm.sendCount_Yz,
			Dm.recvList_x,Dm.recvList_y,Dm.recvList_z,Dm.recvList_X,Dm.recvList_Y,Dm.recvList_Z,
			Dm.recvList_xy,Dm.recvList_XY,Dm.recvList_xY,Dm.recvList_Xy,Dm.recvList_xz,Dm.recvList_XZ,
			Dm.recvList_xZ,Dm.recvList_Xz,Dm.recvList_yz,Dm.recvList_YZ,Dm.recvList_yZ,Dm.recvList_Yz,
			Dm.recvCount_x,Dm.recvCount_y,Dm.recvCount_z,Dm.recvCount_X,Dm.recvCount_Y,Dm.recvCount_Z,
			Dm.recvCount_xy,Dm.recvCount_XY,Dm.recvCount_xY,Dm.recvCount_Xy,Dm.recvCount_xz,Dm.recvCount_XZ,
			Dm.recvCount_xZ,Dm.recvCount_Xz,Dm.recvCount_yz,Dm.recvCount_YZ,Dm.recvCount_yZ,Dm.recvCount_Yz,
			Dm.rank_x,Dm.rank_y,Dm.rank_z,Dm.rank_X,Dm.rank_Y,Dm.rank_Z,Dm.rank_xy,Dm.rank_XY,Dm.rank_xY,
			Dm.rank_Xy,Dm.rank_xz,Dm.rank_XZ,Dm.rank_xZ,Dm.rank_Xz,Dm.rank_yz,Dm.rank_YZ,Dm.rank_yZ,Dm.rank_Yz);
	//...........................................................................
	CommunicateMeshHalo(SDn_z, Dm.Comm,
			sendMeshData_x,sendMeshData_y,sendMeshData_z,sendMeshData_X,sendMeshData_Y,sendMeshData_Z,
			sendMeshData_xy,sendMeshData_XY,sendMeshData_xY,sendMeshData_Xy,sendMeshData_xz,sendMeshData_XZ,
			sendMeshData_xZ,sendMeshData_Xz,sendMeshData_yz,sendMeshData_YZ,sendMeshData_yZ,sendMeshData_Yz,
			recvMeshData_x,recvMeshData_y,recvMeshData_z,recvMeshData_X,recvMeshData_Y,recvMeshData_Z,
			recvMeshData_xy,recvMeshData_XY,recvMeshData_xY,recvMeshData_Xy,recvMeshData_xz,recvMeshData_XZ,
			recvMeshData_xZ,recvMeshData_Xz,recvMeshData_yz,recvMeshData_YZ,recvMeshData_yZ,recvMeshData_Yz,
			Dm.sendList_x,Dm.sendList_y,Dm.sendList_z,Dm.sendList_X,Dm.sendList_Y,Dm.sendList_Z,
			Dm.sendList_xy,Dm.sendList_XY,Dm.sendList_xY,Dm.sendList_Xy,Dm.sendList_xz,Dm.sendList_XZ,
			Dm.sendList_xZ,Dm.sendList_Xz,Dm.sendList_yz,Dm.sendList_YZ,Dm.sendList_yZ,Dm.sendList_Yz,
			Dm.sendCount_x,Dm.sendCount_y,Dm.sendCount_z,Dm.sendCount_X,Dm.sendCount_Y,Dm.sendCount_Z,
			Dm.sendCount_xy,Dm.sendCount_XY,Dm.sendCount_xY,Dm.sendCount_Xy,Dm.sendCount_xz,Dm.sendCount_XZ,
			Dm.sendCount_xZ,Dm.sendCount_Xz,Dm.sendCount_yz,Dm.sendCount_YZ,Dm.sendCount_yZ,Dm.sendCount_Yz,
			Dm.recvList_x,Dm.recvList_y,Dm.recvList_z,Dm.recvList_X,Dm.recvList_Y,Dm.recvList_Z,
			Dm.recvList_xy,Dm.recvList_XY,Dm.recvList_xY,Dm.recvList_Xy,Dm.recvList_xz,Dm.recvList_XZ,
			Dm.recvList_xZ,Dm.recvList_Xz,Dm.recvList_yz,Dm.recvList_YZ,Dm.recvList_yZ,Dm.recvList_Yz,
			Dm.recvCount_x,Dm.recvCount_y,Dm.recvCount_z,Dm.recvCount_X,Dm.recvCount_Y,Dm.recvCount_Z,
			Dm.recvCount_xy,Dm.recvCount_XY,Dm.recvCount_xY,Dm.recvCount_Xy,Dm.recvCount_xz,Dm.recvCount_XZ,
			Dm.recvCount_xZ,Dm.recvCount_Xz,Dm.recvCount_yz,Dm.recvCount_YZ,Dm.recvCount_yZ,Dm.recvCount_Yz,
			Dm.rank_x,Dm.rank_y,Dm.rank_z,Dm.rank_X,Dm.rank_Y,Dm.rank_Z,Dm.rank_xy,Dm.rank_XY,Dm.rank_xY,
			Dm.rank_Xy,Dm.rank_xz,Dm.rank_XZ,Dm.rank_xZ,Dm.rank_Xz,Dm.rank_yz,Dm.rank_YZ,Dm.rank_yZ,Dm.rank_Yz);
	//...........................................................................	//...........................................................................
	CommunicateMeshHalo(DelPhi, Dm.Comm,
			sendMeshData_x,sendMeshData_y,sendMeshData_z,sendMeshData_X,sendMeshData_Y,sendMeshData_Z,
			sendMeshData_xy,sendMeshData_XY,sendMeshData_xY,sendMeshData_Xy,sendMeshData_xz,sendMeshData_XZ,
			sendMeshData_xZ,sendMeshData_Xz,sendMeshData_yz,sendMeshData_YZ,sendMeshData_yZ,sendMeshData_Yz,
			recvMeshData_x,recvMeshData_y,recvMeshData_z,recvMeshData_X,recvMeshData_Y,recvMeshData_Z,
			recvMeshData_xy,recvMeshData_XY,recvMeshData_xY,recvMeshData_Xy,recvMeshData_xz,recvMeshData_XZ,
			recvMeshData_xZ,recvMeshData_Xz,recvMeshData_yz,recvMeshData_YZ,recvMeshData_yZ,recvMeshData_Yz,
			Dm.sendList_x,Dm.sendList_y,Dm.sendList_z,Dm.sendList_X,Dm.sendList_Y,Dm.sendList_Z,
			Dm.sendList_xy,Dm.sendList_XY,Dm.sendList_xY,Dm.sendList_Xy,Dm.sendList_xz,Dm.sendList_XZ,
			Dm.sendList_xZ,Dm.sendList_Xz,Dm.sendList_yz,Dm.sendList_YZ,Dm.sendList_yZ,Dm.sendList_Yz,
			Dm.sendCount_x,Dm.sendCount_y,Dm.sendCount_z,Dm.sendCount_X,Dm.sendCount_Y,Dm.sendCount_Z,
			Dm.sendCount_xy,Dm.sendCount_XY,Dm.sendCount_xY,Dm.sendCount_Xy,Dm.sendCount_xz,Dm.sendCount_XZ,
			Dm.sendCount_xZ,Dm.sendCount_Xz,Dm.sendCount_yz,Dm.sendCount_YZ,Dm.sendCount_yZ,Dm.sendCount_Yz,
			Dm.recvList_x,Dm.recvList_y,Dm.recvList_z,Dm.recvList_X,Dm.recvList_Y,Dm.recvList_Z,
			Dm.recvList_xy,Dm.recvList_XY,Dm.recvList_xY,Dm.recvList_Xy,Dm.recvList_xz,Dm.recvList_XZ,
			Dm.recvList_xZ,Dm.recvList_Xz,Dm.recvList_yz,Dm.recvList_YZ,Dm.recvList_yZ,Dm.recvList_Yz,
			Dm.recvCount_x,Dm.recvCount_y,Dm.recvCount_z,Dm.recvCount_X,Dm.recvCount_Y,Dm.recvCount_Z,
			Dm.recvCount_xy,Dm.recvCount_XY,Dm.recvCount_xY,Dm.recvCount_Xy,Dm.recvCount_xz,Dm.recvCount_XZ,
			Dm.recvCount_xZ,Dm.recvCount_Xz,Dm.recvCount_yz,Dm.recvCount_YZ,Dm.recvCount_yZ,Dm.recvCount_Yz,
			Dm.rank_x,Dm.rank_y,Dm.rank_z,Dm.rank_X,Dm.rank_Y,Dm.rank_Z,Dm.rank_xy,Dm.rank_XY,Dm.rank_xY,
			Dm.rank_Xy,Dm.rank_xz,Dm.rank_XZ,Dm.rank_xZ,Dm.rank_Xz,Dm.rank_yz,Dm.rank_YZ,Dm.rank_yZ,Dm.rank_Yz);
	//...........................................................................

}

void TwoPhase::ComputeLocal(){
	int i,j,k,n;
	double delphi;
	for (int c=0;c<ncubes;c++){
		// Get cube from the list
		i = cubeList(0,c);
		j = cubeList(1,c);
		k = cubeList(2,c);

		n_nw_pts=n_ns_pts=n_ws_pts=n_nws_pts=n_local_sol_pts=n_local_nws_pts=0;
		n_nw_tris=n_ns_tris=n_ws_tris=n_nws_seg=n_local_sol_tris=0;
		//...........................................................................
		//...........................................................................
		// Compute volume averages
		for (int p=0;p<8;p++){

			if ( SDs(i+cube[p][0],j+cube[p][1],k+cube[p][2]) > 0 ){
				// 1-D index for this cube corner
				n = i+cube[p][0] + (j+cube[p][1])*Nx + (k+cube[p][2])*Nx*Ny;
				// compute the norm of the gradient of the phase indicator field
				// Compute the non-wetting phase volume contribution
				if ( Phase(i+cube[p][0],j+cube[p][1],k+cube[p][2]) > 0 ){
					nwp_volume += 0.125;
					// volume the excludes the interfacial region
					if (DelPhi.data[n] < 1e-4){
						vol_n += 0.125;
						// pressure
						pan += 0.125*Press.data[n];
						// velocity
						van(0) += 0.125*Vel_x.data[n];
						van(1) += 0.125*Vel_y.data[n];
						van(2) += 0.125*Vel_z.data[n];
					}
				}
				else{
					wp_volume += 0.125;
					if (DelPhi.data[n] < 1e-4){
						// volume the excludes the interfacial region
						vol_w += 0.125;
						// pressure
						paw += 0.125*Press.data[n];
						// velocity
						vaw(0) += 0.125*Vel_x.data[n];
						vaw(1) += 0.125*Vel_y.data[n];
						vaw(2) += 0.125*Vel_z.data[n];
					}
				}
			}
		}

		//...........................................................................
		// Construct the interfaces and common curve
		pmmc_ConstructLocalCube(SDs, SDn, solid_isovalue, fluid_isovalue,
				nw_pts, nw_tris, Values, ns_pts, ns_tris, ws_pts, ws_tris,
				local_nws_pts, nws_pts, nws_seg, local_sol_pts, local_sol_tris,
				n_local_sol_tris, n_local_sol_pts, n_nw_pts, n_nw_tris,
				n_ws_pts, n_ws_tris, n_ns_tris, n_ns_pts, n_local_nws_pts, n_nws_pts, n_nws_seg,
				i, j, k, Nx, Ny, Nz);

		// Integrate the contact angle
		efawns += pmmc_CubeContactAngle(CubeValues,Values,SDn_x,SDn_y,SDn_z,SDs_x,SDs_y,SDs_z,
				local_nws_pts,i,j,k,n_local_nws_pts);

		// Integrate the mean curvature
		Jwn    += pmmc_CubeSurfaceInterpValue(CubeValues,MeanCurvature,nw_pts,nw_tris,Values,i,j,k,n_nw_pts,n_nw_tris);
		Kwn    += pmmc_CubeSurfaceInterpValue(CubeValues,GaussCurvature,nw_pts,nw_tris,Values,i,j,k,n_nw_pts,n_nw_tris);

		// Integrate the trimmed mean curvature (hard-coded to use a distance of 4 pixels)
		pmmc_CubeTrimSurfaceInterpValues(CubeValues,MeanCurvature,SDs,nw_pts,nw_tris,Values,DistanceValues,
				i,j,k,n_nw_pts,n_nw_tris,trimdist,trawn,trJwn);

		pmmc_CubeTrimSurfaceInterpInverseValues(CubeValues,MeanCurvature,SDs,nw_pts,nw_tris,Values,DistanceValues,
				i,j,k,n_nw_pts,n_nw_tris,trimdist,dummy,trRwn);

		// Compute the normal speed of the interface
		pmmc_InterfaceSpeed(dPdt, SDn_x, SDn_y, SDn_z, CubeValues, nw_pts, nw_tris,
				NormalVector, InterfaceSpeed, vawn, i, j, k, n_nw_pts, n_nw_tris);

		pmmc_CommonCurveSpeed(CubeValues, dPdt, vawns, SDn_x, SDn_y, SDn_z,SDs_x,SDs_y,SDs_z,
				local_nws_pts,i,j,k,n_local_nws_pts);

		pmmc_CurveCurvature(SDn, SDs, KNwns_values, KGwns_values, KNwns, KGwns,
				nws_pts, n_nws_pts, i, j, k);

		As  += pmmc_CubeSurfaceArea(local_sol_pts,local_sol_tris,n_local_sol_tris);

		// Compute the surface orientation and the interfacial area
		awn += pmmc_CubeSurfaceOrientation(Gwn,nw_pts,nw_tris,n_nw_tris);
		ans += pmmc_CubeSurfaceOrientation(Gns,ns_pts,ns_tris,n_ns_tris);
		aws += pmmc_CubeSurfaceOrientation(Gws,ws_pts,ws_tris,n_ws_tris);
		lwns +=  pmmc_CubeCurveLength(local_nws_pts,n_local_nws_pts);
		//...........................................................................
	}
}

void TwoPhase::Reduce(){
	int i;
	double iVol_global=1.0/Volume;
	//...........................................................................
	MPI_Barrier(Dm.Comm);
	MPI_Allreduce(&nwp_volume,&nwp_volume_global,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
	MPI_Allreduce(&wp_volume,&wp_volume_global,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
	MPI_Allreduce(&awn,&awn_global,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
	MPI_Allreduce(&ans,&ans_global,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
	MPI_Allreduce(&aws,&aws_global,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
	MPI_Allreduce(&lwns,&lwns_global,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
	MPI_Allreduce(&As,&As_global,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
	MPI_Allreduce(&Jwn,&Jwn_global,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
	MPI_Allreduce(&Kwn,&Kwn_global,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
	MPI_Allreduce(&KGwns,&KGwns_global,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
	MPI_Allreduce(&KNwns,&KNwns_global,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
	MPI_Allreduce(&efawns,&efawns_global,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
	// Phase averages
	MPI_Allreduce(&vol_w,&vol_w_global,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
	MPI_Allreduce(&vol_n,&vol_n_global,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
	MPI_Allreduce(&paw,&paw_global,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
	MPI_Allreduce(&pan,&pan_global,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
	MPI_Allreduce(&vaw(0),&vaw_global(0),3,MPI_DOUBLE,MPI_SUM,Dm.Comm);
	MPI_Allreduce(&van(0),&van_global(0),3,MPI_DOUBLE,MPI_SUM,Dm.Comm);
	MPI_Allreduce(&vawn(0),&vawn_global(0),3,MPI_DOUBLE,MPI_SUM,Dm.Comm);
	MPI_Allreduce(&vawns(0),&vawns_global(0),3,MPI_DOUBLE,MPI_SUM,Dm.Comm);
	MPI_Allreduce(&Gwn(0),&Gwn_global(0),6,MPI_DOUBLE,MPI_SUM,Dm.Comm);
	MPI_Allreduce(&Gns(0),&Gns_global(0),6,MPI_DOUBLE,MPI_SUM,Dm.Comm);
	MPI_Allreduce(&Gws(0),&Gws_global(0),6,MPI_DOUBLE,MPI_SUM,Dm.Comm);
	MPI_Allreduce(&trawn,&trawn_global,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
	MPI_Allreduce(&trJwn,&trJwn_global,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
	MPI_Allreduce(&trRwn,&trRwn_global,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
	MPI_Barrier(Dm.Comm);

	// Normalize the phase averages
	// (density of both components = 1.0)
	if (vol_w_global > 0.0){
		paw_global = paw_global / vol_w_global;
		vaw_global(0) = vaw_global(0) / vol_w_global;
		vaw_global(1) = vaw_global(1) / vol_w_global;
		vaw_global(2) = vaw_global(2) / vol_w_global;
	}
	if (vol_n_global > 0.0){
		pan_global = pan_global / vol_n_global;
		van_global(0) = van_global(0) / vol_n_global;
		van_global(1) = van_global(1) / vol_n_global;
		van_global(2) = van_global(2) / vol_n_global;
	}
	// Normalize surface averages by the interfacial area
	if (awn_global > 0.0){
		Jwn_global /= awn_global;
		Kwn_global /= awn_global;
		for (i=0; i<3; i++) vawn_global(i) /= awn_global;
		for (i=0; i<6; i++)	Gwn_global(i) /= awn_global;
	}
	if (lwns_global > 0.0){
		efawns_global /= lwns_global;
		KNwns_global /= lwns_global;
		KGwns_global /= lwns_global;
		for (i=0; i<3; i++)	vawns_global(i) /= lwns_global;
	}
	if (trawn_global > 0.0){
		trJwn_global /= trawn_global;
		trRwn_global /= trawn_global;
		trRwn_global = 2.0*fabs(trRwn_global);
		trJwn_global = fabs(trJwn_global);
	}

	if (ans_global > 0.0)	for (i=0; i<6; i++)		Gns_global(i) /= ans_global;
	if (aws_global > 0.0)	for (i=0; i<6; i++)		Gws_global(i) /= aws_global;

	//sat_w = 1.0 - nwp_volume_global*iVol_global/porosity;
	sat_w = 1.0 - nwp_volume_global/(nwp_volume_global+wp_volume_global);
	// Compute the specific interfacial areas and common line length (dimensionless per unit volume)
	awn_global = awn_global*iVol_global;
	ans_global = ans_global*iVol_global;
	aws_global = aws_global*iVol_global;
	dEs = dEs*iVol_global;
	lwns_global = lwns_global*iVol_global;
}

void TwoPhase::NonDimensionalize(double D, double viscosity, double IFT){
	awn_global *= D;
	ans_global *= D;
	ans_global *= D;
	lwns_global *= D*D;
}

void TwoPhase::PrintAll(int timestep){
	if (Dm.rank==0){
		fprintf(TIMELOG,"%i %.5g ",timestep-5,dEs);										// change in surface energy
		fprintf(TIMELOG,"%.5g %.5g %.5g ",sat_w,paw_global,pan_global);					// saturation and pressure
		fprintf(TIMELOG,"%.5g %.5g %.5g ",awn_global,ans_global,aws_global);				// interfacial areas
		fprintf(TIMELOG,"%.5g %.5g ",Jwn_global, Kwn_global);								// curvature of wn interface
		fprintf(TIMELOG,"%.5g ",lwns_global);											// common curve length
		fprintf(TIMELOG,"%.5g ",efawns_global);											// average contact angle
		fprintf(TIMELOG,"%.5g %.5g ",KNwns_global, KGwns_global);								// curvature of wn interface
		fprintf(TIMELOG,"%.5g %.5g %.5g ",vaw_global(0),vaw_global(1),vaw_global(2));	// average velocity of w phase
		fprintf(TIMELOG,"%.5g %.5g %.5g ",van_global(0),van_global(1),van_global(2));	// average velocity of n phase
		fprintf(TIMELOG,"%.5g %.5g %.5g ",vawn_global(0),vawn_global(1),vawn_global(2));	// velocity of wn interface
		fprintf(TIMELOG,"%.5g %.5g %.5g ",vawns_global(0),vawns_global(1),vawns_global(2));	// velocity of wn interface
		fprintf(TIMELOG,"%.5g %.5g %.5g %.5g %.5g %.5g ",
				Gwn_global(0),Gwn_global(1),Gwn_global(2),Gwn_global(3),Gwn_global(4),Gwn_global(5));	// orientation of wn interface
		fprintf(TIMELOG,"%.5g %.5g %.5g %.5g %.5g %.5g ",
				Gns_global(0),Gns_global(1),Gns_global(2),Gns_global(3),Gns_global(4),Gns_global(5));	// orientation of ns interface
		fprintf(TIMELOG,"%.5g %.5g %.5g %.5g %.5g %.5g ",
				Gws_global(0),Gws_global(1),Gws_global(2),Gws_global(3),Gws_global(4),Gws_global(5));	// orientation of ws interface
		fprintf(TIMELOG,"%.5g %.5g %.5g\n",trawn_global, trJwn_global, trRwn_global);						// Trimmed curvature
		fflush(TIMELOG);
	}
}


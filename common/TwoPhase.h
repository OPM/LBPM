// Header file for two-phase averaging class
#include <vector>

#include "pmmc.h"
#include "Domain.h"
#include "Communication.h"
#include "analysis/analysis.h"

#include "shared_ptr.h"
#include "common/Utilities.h"
#include "common/MPI_Helpers.h"
#include "IO/MeshDatabase.h"
#include "IO/Reader.h"
#include "IO/Writer.h"

#define BLOB_AVG_COUNT 33

// Array access for averages defined by the following
#define VOL 0
#define TRIMVOL 1
#define PRS 2
#define AWN 3
#define AWS 4
#define ANS 5
#define LWNS 6
#define JWN 7
#define KWN 8
#define CWNS 9
#define KNWNS 10
#define KGWNS 11
#define VX 12
#define VY 13
#define VZ 14
#define VSQ 15
#define VWNX 16
#define VWNY 17
#define VWNZ 18
#define VWNSX 19
#define VWNSY 20
#define VWNSZ 21
#define GWNXX 22
#define GWNYY 23
#define GWNZZ 24
#define GWNXY 25
#define GWNXZ 26
#define GWNYZ 27
#define TRAWN 28
#define TRJWN 29
#define CMX 30
#define CMY 31
#define CMZ 32

class TwoPhase{

	//...........................................................................
	int n_nw_pts,n_ns_pts,n_ws_pts,n_nws_pts,n_local_sol_pts,n_local_nws_pts;
	int n_nw_tris,n_ns_tris,n_ws_tris,n_nws_seg,n_local_sol_tris;
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

	// Temporary storage arrays
	DoubleArray CubeValues;
	DoubleArray Values;
	DoubleArray DistanceValues;
	DoubleArray KGwns_values;
	DoubleArray KNwns_values;
	DoubleArray InterfaceSpeed;
	DoubleArray NormalVector;

	DoubleArray RecvBuffer;

	// CSV / text file where time history of averages is saved
	FILE *TIMELOG;
	FILE *NWPLOG;
	FILE *WPLOG;

public:
	//...........................................................................
	Domain& Dm;
	int NumberComponents_WP,NumberComponents_NWP;
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
	IntArray PhaseID;	// Phase ID array (solid=0, non-wetting=1, wetting=2)
	IntArray Label_WP;
	IntArray Label_NWP;
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
	//	Container for averages;
	DoubleArray ComponentAverages_WP;
	DoubleArray ComponentAverages_NWP;
	//...........................................................................
	TwoPhase(Domain &dm) : Dm(dm){
		Nx=dm.Nx; Ny=dm.Ny; Nz=dm.Nz;
		Volume=(Nx-2)*(Ny-2)*(Nz-2)*Dm.nprocx*Dm.nprocy*Dm.nprocz*1.0;

		// Global arrays
		PhaseID.resize(Nx,Ny,Nz);
		Label_WP.resize(Nx,Ny,Nz);
		Label_NWP.resize(Nx,Ny,Nz);
		SDn.resize(Nx,Ny,Nz);
		SDs.resize(Nx,Ny,Nz);
		Phase.resize(Nx,Ny,Nz);
		Press.resize(Nx,Ny,Nz);
		dPdt.resize(Nx,Ny,Nz);
		MeanCurvature.resize(Nx,Ny,Nz);
		GaussCurvature.resize(Nx,Ny,Nz);
		SDs_x.resize(Nx,Ny,Nz);		// Gradient of the signed distance
		SDs_y.resize(Nx,Ny,Nz);
		SDs_z.resize(Nx,Ny,Nz);
		SDn_x.resize(Nx,Ny,Nz);		// Gradient of the signed distance
		SDn_y.resize(Nx,Ny,Nz);
		SDn_z.resize(Nx,Ny,Nz);
		DelPhi.resize(Nx,Ny,Nz);
		Phase_tplus.resize(Nx,Ny,Nz);
		Phase_tminus.resize(Nx,Ny,Nz);
		Vel_x.resize(Nx,Ny,Nz);		// Gradient of the phase indicator field
		Vel_y.resize(Nx,Ny,Nz);
		Vel_z.resize(Nx,Ny,Nz);
		//.........................................
		// Allocate cube storage space
		CubeValues.resize(2,2,2);
		nw_tris.resize(3,20);
		ns_tris.resize(3,20);
		ws_tris.resize(3,20);
		nws_seg.resize(2,20);
		local_sol_tris.resize(3,18);
		nw_pts=DTMutableList<Point>(20);
		ns_pts=DTMutableList<Point>(20);
		ws_pts=DTMutableList<Point>(20);
		nws_pts=DTMutableList<Point>(20);
		local_nws_pts=DTMutableList<Point>(20);
		local_sol_pts=DTMutableList<Point>(20);
		tmp=DTMutableList<Point>(20);
		//.........................................
		Values.resize(20);
		DistanceValues.resize(20);
		KGwns_values.resize(20);
		KNwns_values.resize(20);
		InterfaceSpeed.resize(20);
		NormalVector.resize(60);
		//.........................................
		van.resize(3);
		vaw.resize(3);
		vawn.resize(3);
		vawns.resize(3);
		Gwn.resize(6);
		Gns.resize(6);
		Gws.resize(6);
		van_global.resize(3);
		vaw_global.resize(3);
		vawn_global.resize(3);
		vawns_global.resize(3);
		Gwn_global.resize(6);
		Gns_global.resize(6);
		Gws_global.resize(6);
		//.........................................
		if (Dm.rank==0){
			TIMELOG = fopen("timelog.tcat","a+");
			if (fseek(TIMELOG,0,SEEK_SET) == fseek(TIMELOG,0,SEEK_CUR)){
				// If timelog is empty, write a short header to list the averages
				//fprintf(TIMELOG,"--------------------------------------------------------------------------------------\n");
				fprintf(TIMELOG,"time dEs ");								// Timestep, Change in Surface Energy
				fprintf(TIMELOG,"sw pw pn awn ans aws Jwn Kwn lwns cwns KNwns KGwns ");	// Scalar averages
				fprintf(TIMELOG,"vawx vawy vawz vanx vany vanz ");			// Velocity averages
				fprintf(TIMELOG,"vawnx vawny vawnz vawnsx vawnsy vawnsz ");
				fprintf(TIMELOG,"Gwnxx Gwnyy Gwnzz Gwnxy Gwnxz Gwnyz ");				// Orientation tensors
				fprintf(TIMELOG,"Gwsxx Gwsyy Gwszz Gwsxy Gwsxz Gwsyz ");
				fprintf(TIMELOG,"Gnsxx Gnsyy Gnszz Gnsxy Gnsxz Gnsyz ");
				fprintf(TIMELOG,"trawn trJwn trRwn\n");								// trimmed curvature for wn surface
				//fprintf(TIMELOG,"--------------------------------------------------------------------------------------\n");
			}

			NWPLOG = fopen("components.NWP.tcat","a+");
			fprintf(NWPLOG,"time label vol pn awn ans Jwn Kwn lwns cwns ");
			fprintf(NWPLOG,"vx vy vz vwnx vwny vwnz vwnsx vwnsy vwnsz vsq ");
			fprintf(NWPLOG,"Gwnxx Gwnyy Gwnzz Gwnxy Gwnxz Gwnyz trawn trJwn\n");

			WPLOG = fopen("components.WP.tcat","a+");
			fprintf(WPLOG,"time label vol pw awn ans Jwn Kwn lwns cwns ");
			fprintf(WPLOG,"vx vy vz vwnx vwny vwnz vwnsx vwnsy vwnsz vsq ");
			fprintf(WPLOG,"Gwnxx Gwnyy Gwnzz Gwnxy Gwnxz Gwnyz trawn trJwn\n");
		}
	}
	~TwoPhase(){

	}

	void Initialize();
//	void SetupCubes(Domain &Dm);
	void UpdateMeshValues();
	void UpdateSolid();
	void ComputeDelPhi();
	void ColorToSignedDistance(double Beta, DoubleArray &ColorData, DoubleArray &DistData);
	void ComputeLocal();
	void ComponentAverages();
	void Reduce();
	void WriteSurfaces(int logcount);
	void NonDimensionalize(double D, double viscosity, double IFT);
	void PrintAll(int timestep);
	int GetCubeLabel(int i, int j, int k, IntArray &BlobLabel);
	void SortBlobs();
	void PrintComponents(int timestep);
};

void TwoPhase::ColorToSignedDistance(double Beta, DoubleArray &ColorData, DoubleArray &DistData){

	double factor,temp,value;
	factor=0.5/Beta;
/*	for (int n=0; n<Nx*Ny*Nz; n++){
	 	value = ColorData[n];
		if (value > 0.999 ) DistData[n] = 4.0;
		else if (value < -0.999 ) DistData[n] = -4.0;
		else 	DistData[n] = factor*log((1.0+value)/(1.0-value));
		if (DistData[n] > 1.0)  DistData[n] = 1.0;
		if (DistData[n] < -1.0) DistData[n] = -1.0;
	}
	// Initialize to -1,1 (segmentation)
	for (int k=0; k<Nz; k++){
		for (int j=0; j<Ny; j++){
			for (int i=0; i<Nx; i++){

				value = ColorData(i,j,k);
				temp = factor*log((1.0+value)/(1.0-value));
				if (temp > 1.0) DistData(i,j,k) = 1.0;
				else if (temp < -1.0) DistData(i,j,k) = -1.0;
				else DistData(i,j,k) = temp;
			}
		}
	}

	SSO(DistData,Dm.id,Dm,10);
*/
	for (int k=0; k<Nz; k++){
		for (int j=0; j<Ny; j++){
			for (int i=0; i<Nx; i++){
				DistData(i,j,k) = ColorData(i,j,k);
			}
		}
	}
//	for (int n=0; n<Nx*Ny*Nz; n++)	DistData[n] = ColorData[n];
}

void TwoPhase::ComputeDelPhi(){

	int i,j,k;
	double fx,fy,fz;

	Dm.CommunicateMeshHalo(Phase);
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
/*
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
*/
void TwoPhase::UpdateSolid(){
	Dm.CommunicateMeshHalo(SDs);
	//...........................................................................
	// Gradient of the Signed Distance function
	//...........................................................................
	pmmc_MeshGradient(SDs,SDs_x,SDs_y,SDs_z,Nx,Ny,Nz);
	//...........................................................................
	Dm.CommunicateMeshHalo(SDs_x);
	//...........................................................................
	Dm.CommunicateMeshHalo(SDs_y);
	//...........................................................................
	Dm.CommunicateMeshHalo(SDs_z);
	//...........................................................................
}

void TwoPhase::UpdateMeshValues(){
	int i,j,k,n;
	//...........................................................................
	// Compute the gradients of the phase indicator and signed distance fields
	pmmc_MeshGradient(SDn,SDn_x,SDn_y,SDn_z,Nx,Ny,Nz);
	//...........................................................................
	// Gradient of the phase indicator field
	//...........................................................................
	Dm.CommunicateMeshHalo(SDn_x);
	//...........................................................................
	Dm.CommunicateMeshHalo(SDn_y);
	//...........................................................................
	Dm.CommunicateMeshHalo(SDn_z);
	//...........................................................................
	// Compute the mesh curvature of the phase indicator field
	pmmc_MeshCurvature(SDn, MeanCurvature, GaussCurvature, Nx, Ny, Nz);
	//...........................................................................
	// Update the time derivative of non-dimensional density field
	// Map Phase_tplus and Phase_tminus
	for (int n=0; n<Nx*Ny*Nz; n++)	dPdt(n) = 0.1*(Phase_tplus(n) - Phase_tminus(n));
	//...........................................................................
	Dm.CommunicateMeshHalo(Press);
	//...........................................................................
	Dm.CommunicateMeshHalo(Vel_x);
	//...........................................................................
	Dm.CommunicateMeshHalo(Vel_y);
	//...........................................................................
	Dm.CommunicateMeshHalo(Vel_z);
	//...........................................................................
	Dm.CommunicateMeshHalo(MeanCurvature);
	//...........................................................................
	Dm.CommunicateMeshHalo(GaussCurvature);
	//...........................................................................
	Dm.CommunicateMeshHalo(DelPhi);
	//...........................................................................
	// Initializing the blob ID
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				n = k*Nx*Ny+j*Nx+i;
				if (Dm.id[n] == 0){
					// Solid phase
					PhaseID(i,j,k) = 0;
				}
				else if (SDn(i,j,k) < 0.0){
					// wetting phase
					PhaseID(i,j,k) = 2;
				}
				else {
					// non-wetting phase
					PhaseID(i,j,k) = 1;
				}
			}
		}
	}

}
void TwoPhase::ComputeLocal(){
	int i,j,k,n,kmin,kmax;
	int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};

	// If external boundary conditions are set, do not average over the inlet
	kmin=1; kmax=Nz-1;
	if (Dm.BoundaryCondition > 0 && Dm.kproc == 0) kmin=4;
	if (Dm.BoundaryCondition > 0 && Dm.kproc == Dm.nprocz-1) kmax=Nz-4;

	for (k=kmin; k<kmax; k++){
		for (j=1; j<Ny-1; j++){
			for (i=1; i<Nx-1; i++){
				//...........................................................................
				n_nw_pts=n_ns_pts=n_ws_pts=n_nws_pts=n_local_sol_pts=n_local_nws_pts=0;
				n_nw_tris=n_ns_tris=n_ws_tris=n_nws_seg=n_local_sol_tris=0;
				//...........................................................................
				// Compute volume averages
				for (int p=0;p<8;p++){
					n = i+cube[p][0] + (j+cube[p][1])*Nx + (k+cube[p][2])*Nx*Ny;
					if ( Dm.id[n] != 0 ){
						// 1-D index for this cube corner
						// compute the norm of the gradient of the phase indicator field
						// Compute the non-wetting phase volume contribution
						if ( Phase(i+cube[p][0],j+cube[p][1],k+cube[p][2]) > 0 ){
							nwp_volume += 0.125;
							// velocity
							van(0) += 0.125*Vel_x(n);
							van(1) += 0.125*Vel_y(n);
							van(2) += 0.125*Vel_z(n);
							// volume the excludes the interfacial region
							if (DelPhi(n) < 1e-4){
								vol_n += 0.125;
								// pressure
								pan += 0.125*Press(n);

							}
						}
						else{
							wp_volume += 0.125;
							// velocity
							vaw(0) += 0.125*Vel_x(n);
							vaw(1) += 0.125*Vel_y(n);
							vaw(2) += 0.125*Vel_z(n);
							if (DelPhi(n) < 1e-4){
								// volume the excludes the interfacial region
								vol_w += 0.125;
								// pressure
								paw += 0.125*Press(n);

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

				// wn interface averages
				if (n_nw_pts > 0){
					awn += pmmc_CubeSurfaceOrientation(Gwn,nw_pts,nw_tris,n_nw_tris);
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
				}
				// wns common curve averages
				if (n_local_nws_pts > 0){
					efawns += pmmc_CubeContactAngle(CubeValues,Values,SDn_x,SDn_y,SDn_z,SDs_x,SDs_y,SDs_z,
							local_nws_pts,i,j,k,n_local_nws_pts);

					pmmc_CommonCurveSpeed(CubeValues, dPdt, vawns, SDn_x, SDn_y, SDn_z,SDs_x,SDs_y,SDs_z,
							local_nws_pts,i,j,k,n_local_nws_pts);

					pmmc_CurveCurvature(SDn, SDs, SDn_x, SDn_y, SDn_z, SDs_x, SDs_y,
							SDs_z, KNwns_values, KGwns_values, KNwns, KGwns,
							nws_pts, n_nws_pts, i, j, k);

					lwns +=  pmmc_CubeCurveLength(local_nws_pts,n_local_nws_pts);
				}

				// Solid interface averagees
				if (n_local_sol_tris > 0){
					As  += pmmc_CubeSurfaceArea(local_sol_pts,local_sol_tris,n_local_sol_tris);

					// Compute the surface orientation and the interfacial area
					ans += pmmc_CubeSurfaceOrientation(Gns,ns_pts,ns_tris,n_ns_tris);
					aws += pmmc_CubeSurfaceOrientation(Gws,ws_pts,ws_tris,n_ws_tris);
				}
				//...........................................................................
			}
		}
	}
}

void TwoPhase::ComponentAverages(){
    int i,j,k,n;
    int kmin,kmax;
	int LabelWP,LabelNWP;
	double TempLocal;

	int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};

	LabelNWP=1; LabelWP=2;
	// NOTE: labeling the wetting phase components is tricky! One sandstone media had over 800,000 components
	//NumberComponents_WP = ComputeGlobalPhaseComponent(Dm.Nx-2,Dm.Ny-2,Dm.Nz-2,Dm.rank_info,PhaseID,LabelWP,Label_WP);
	// treat all wetting phase is connected
	NumberComponents_WP=1;
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				Label_WP(i,j,k) = 0;
			}
		}
	}

	// Fewer non-wetting phase features are present
	NumberComponents_NWP = ComputeGlobalPhaseComponent(Dm.Nx-2,Dm.Ny-2,Dm.Nz-2,Dm.rank_info,PhaseID,LabelNWP,Label_NWP);
	//NumberComponents_NWP = ComputeGlobalBlobIDs(Dm.Nx-2,Dm.Ny-2,Dm.Nz-2,Dm.rank_info,SDs,SDn,solid_isovalue,fluid_isovalue,Label_NWP);

	ComponentAverages_WP.resize(BLOB_AVG_COUNT,NumberComponents_WP);
	ComponentAverages_NWP.resize(BLOB_AVG_COUNT,NumberComponents_NWP);

	ComponentAverages_WP.fill(0.0);
	ComponentAverages_NWP.fill(0.0);
	
	if (Dm.rank==0){
		printf("Number of wetting phase components is %i \n",NumberComponents_WP);
		printf("Number of non-wetting phase components is %i \n",NumberComponents_NWP);
	}

	// If external boundary conditions are set, do not average over the inlet
	kmin=1; kmax=Nz-1;
	if (Dm.BoundaryCondition > 0 && Dm.kproc == 0) kmin=4;
	if (Dm.BoundaryCondition > 0 && Dm.kproc == Dm.nprocz-1) kmax=Nz-4;

	for (k=kmin; k<kmax; k++){
		for (j=1; j<Ny-1; j++){
			for (i=1; i<Nx-1; i++){

				LabelWP=GetCubeLabel(i,j,k,Label_WP);
				LabelNWP=GetCubeLabel(i,j,k,Label_NWP);
				
				n_nw_pts=n_ns_pts=n_ws_pts=n_nws_pts=n_local_sol_pts=n_local_nws_pts=0;
				n_nw_tris=n_ns_tris=n_ws_tris=n_nws_seg=n_local_sol_tris=0;

				// Initialize the averaged quantities
				awn = aws = ans = lwns = 0.0;
				vawn(0) = vawn(1) = vawn(2) = 0.0;
				vawns(0) = vawns(1) = vawns(2) = 0.0;
				Gwn(0) = Gwn(1) = Gwn(2) = 0.0;
				Gwn(3) = Gwn(4) = Gwn(5) = 0.0;
				Gws(0) = Gws(1) = Gws(2) = 0.0;
				Gws(3) = Gws(4) = Gws(5) = 0.0;
				Gns(0) = Gns(1) = Gns(2) = 0.0;
				Gns(3) = Gns(4) = Gns(5) = 0.0;
				KGwns = KNwns = 0.0;
				Jwn = Kwn = efawns = 0.0;
				trawn=trJwn=0.0;
				//...........................................................................
				//...........................................................................
				// Compute volume averages
				for (int p=0;p<8;p++){
					n = i+cube[p][0] + (j+cube[p][1])*Nx + (k+cube[p][2])*Nx*Ny;
					if ( Dm.id[n] != 0 ){
						// 1-D index for this cube corner
						// compute the norm of the gradient of the phase indicator field
						// Compute the non-wetting phase volume contribution
						if ( Phase(i+cube[p][0],j+cube[p][1],k+cube[p][2]) > 0.0 && !(LabelNWP < 0) ){
							// volume
							ComponentAverages_NWP(VOL,LabelNWP) += 0.125;
							// velocity
							ComponentAverages_NWP(VX,LabelNWP) += 0.125*Vel_x(n);
							ComponentAverages_NWP(VY,LabelNWP) += 0.125*Vel_y(n);
							ComponentAverages_NWP(VZ,LabelNWP) += 0.125*Vel_z(n);
							// center of mass
							ComponentAverages_NWP(CMX,LabelNWP) += 0.125*(i+cube[p][0]);
							ComponentAverages_NWP(CMY,LabelNWP) += 0.125*(j+cube[p][1]);
							ComponentAverages_NWP(CMZ,LabelNWP) += 0.125*(k+cube[p][2]);

							// twice the kinetic energy
							ComponentAverages_NWP(VSQ,LabelNWP) += 0.125*(Vel_x(n)*Vel_x(n)+Vel_y(n)*Vel_y(n)+Vel_z(n)*Vel_z(n));

							// volume the for pressure averaging excludes the interfacial region
							if (DelPhi(n) < 1e-4 ){
								ComponentAverages_NWP(TRIMVOL,LabelNWP) += 0.125;
								ComponentAverages_NWP(PRS,LabelNWP) += 0.125*Press(n);
							}
						}
						else if (!(LabelWP < 0)){
							ComponentAverages_WP(VOL,LabelWP) += 0.125;
							// velocity
							ComponentAverages_WP(VX,LabelWP) += 0.125*Vel_x(n);
							ComponentAverages_WP(VY,LabelWP)+= 0.125*Vel_y(n);
							ComponentAverages_WP(VZ,LabelWP) += 0.125*Vel_z(n);
							// Center of mass
							ComponentAverages_WP(CMX,LabelWP) += 0.125*(i+cube[p][0]);
							ComponentAverages_WP(CMY,LabelWP) += 0.125*(j+cube[p][1]);
							ComponentAverages_WP(CMZ,LabelWP) += 0.125*(k+cube[p][2]);
							// twice the kinetic energy
							ComponentAverages_WP(VSQ,LabelWP) += 0.125*(Vel_x(n)*Vel_x(n)+Vel_y(n)*Vel_y(n)+Vel_z(n)*Vel_z(n));

							// volume the for pressure averaging excludes the interfacial region
							if (DelPhi(n) < 1e-4){
								ComponentAverages_WP(TRIMVOL,LabelWP) += 0.125;
								ComponentAverages_WP(PRS,LabelWP) += 0.125*Press(n);
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

				//...........................................................................
				// wn interface averages
				if (n_nw_pts>0  && LabelNWP >=0 && LabelWP >=0 ){
					// Mean curvature
					TempLocal = pmmc_CubeSurfaceInterpValue(CubeValues,MeanCurvature,nw_pts,nw_tris,Values,i,j,k,n_nw_pts,n_nw_tris);
					ComponentAverages_WP(JWN,LabelWP) += TempLocal;
					ComponentAverages_NWP(JWN,LabelNWP) += TempLocal;

					// Trimmed Mean curvature
					pmmc_CubeTrimSurfaceInterpValues(CubeValues,MeanCurvature,SDs,nw_pts,nw_tris,Values,DistanceValues,
							i,j,k,n_nw_pts,n_nw_tris,trimdist,trawn,trJwn);
					ComponentAverages_WP(TRAWN,LabelWP) += trawn;
					ComponentAverages_WP(TRJWN,LabelWP) += trJwn;
					ComponentAverages_NWP(TRAWN,LabelNWP) += trawn;
					ComponentAverages_NWP(TRJWN,LabelNWP) += trJwn;

					// Gaussian curvature
					TempLocal = pmmc_CubeSurfaceInterpValue(CubeValues,GaussCurvature,nw_pts,nw_tris,Values,i,j,k,n_nw_pts,n_nw_tris);
					ComponentAverages_WP(KWN,LabelWP) += TempLocal;
					ComponentAverages_NWP(KWN,LabelNWP) += TempLocal;

					// Compute the normal speed of the interface
					pmmc_InterfaceSpeed(dPdt, SDn_x, SDn_y, SDn_z, CubeValues, nw_pts, nw_tris,
							NormalVector, InterfaceSpeed, vawn, i, j, k, n_nw_pts, n_nw_tris);
					ComponentAverages_WP(VWNX,LabelWP) += vawn(0);
					ComponentAverages_WP(VWNY,LabelWP) += vawn(1);
					ComponentAverages_WP(VWNZ,LabelWP) += vawn(2);
					ComponentAverages_NWP(VWNX,LabelNWP) += vawn(0);
					ComponentAverages_NWP(VWNY,LabelNWP) += vawn(1);
					ComponentAverages_NWP(VWNZ,LabelNWP) += vawn(2);

					// Interfacial Area
					TempLocal = pmmc_CubeSurfaceOrientation(Gwn,nw_pts,nw_tris,n_nw_tris);
					ComponentAverages_WP(AWN,LabelWP) += TempLocal;
					ComponentAverages_NWP(AWN,LabelNWP) += TempLocal;

					ComponentAverages_WP(GWNXX,LabelWP) += Gwn(0);
					ComponentAverages_WP(GWNYY,LabelWP) += Gwn(1);
					ComponentAverages_WP(GWNZZ,LabelWP) += Gwn(2);
					ComponentAverages_WP(GWNXY,LabelWP) += Gwn(3);
					ComponentAverages_WP(GWNXZ,LabelWP) += Gwn(4);
					ComponentAverages_WP(GWNYZ,LabelWP) += Gwn(5);

					ComponentAverages_NWP(GWNXX,LabelNWP) += Gwn(0);
					ComponentAverages_NWP(GWNYY,LabelNWP) += Gwn(1);
					ComponentAverages_NWP(GWNZZ,LabelNWP) += Gwn(2);
					ComponentAverages_NWP(GWNXY,LabelNWP) += Gwn(3);
					ComponentAverages_NWP(GWNXZ,LabelNWP) += Gwn(4);
					ComponentAverages_NWP(GWNYZ,LabelNWP) += Gwn(5);

				}
				//...........................................................................
				// Common curve averages
				if (n_local_nws_pts > 0  && LabelNWP >=0 && LabelWP >=0){
					// Contact angle
					TempLocal = pmmc_CubeContactAngle(CubeValues,Values,SDn_x,SDn_y,SDn_z,SDs_x,SDs_y,SDs_z,
							local_nws_pts,i,j,k,n_local_nws_pts);
					ComponentAverages_WP(CWNS,LabelWP) += TempLocal;
					ComponentAverages_NWP(CWNS,LabelNWP) += TempLocal;

					// Kinematic velocity of the common curve
					pmmc_CommonCurveSpeed(CubeValues, dPdt, vawns, SDn_x, SDn_y, SDn_z,SDs_x,SDs_y,SDs_z,
							local_nws_pts,i,j,k,n_local_nws_pts);
					ComponentAverages_WP(VWNSX,LabelWP) += vawns(0);
					ComponentAverages_WP(VWNSY,LabelWP) += vawns(1);
					ComponentAverages_WP(VWNSZ,LabelWP) += vawns(2);
					ComponentAverages_NWP(VWNSX,LabelNWP) += vawns(0);
					ComponentAverages_NWP(VWNSY,LabelNWP) += vawns(1);
					ComponentAverages_NWP(VWNSZ,LabelNWP) += vawns(2);

					// Curvature of the common curve
					pmmc_CurveCurvature(SDn, SDs, SDn_x, SDn_y, SDn_z, SDs_x, SDs_y,
							SDs_z, KNwns_values, KGwns_values, KNwns, KGwns,
							nws_pts, n_nws_pts, i, j, k);
					ComponentAverages_WP(KNWNS,LabelWP) += KNwns;
					ComponentAverages_WP(KGWNS,LabelWP) += KGwns;
					ComponentAverages_NWP(KNWNS,LabelNWP) += KNwns;
					ComponentAverages_NWP(KGWNS,LabelNWP) += KGwns;

					// Length of the common curve
					TempLocal = pmmc_CubeCurveLength(local_nws_pts,n_local_nws_pts);
					ComponentAverages_NWP(LWNS,LabelNWP) += TempLocal;
				}
				//...........................................................................
				// Solid interface averages
				if (n_local_sol_pts > 0  && LabelWP >=0){
					As  += pmmc_CubeSurfaceArea(local_sol_pts,local_sol_tris,n_local_sol_tris);
					// Compute the surface orientation and the interfacial area

					TempLocal = pmmc_CubeSurfaceOrientation(Gws,ws_pts,ws_tris,n_ws_tris);
					ComponentAverages_WP(AWS,LabelWP) += TempLocal;
				}
				if (n_ns_pts > 0  && LabelNWP >=0 ){
					TempLocal = pmmc_CubeSurfaceOrientation(Gns,ns_pts,ns_tris,n_ns_tris);
					ComponentAverages_NWP(ANS,LabelNWP) += TempLocal;
				}
				//...........................................................................

			}
		}
	}

	if (Dm.rank==0){
		printf("Component averages computed locally -- reducing result... \n");
	}
	// Globally reduce the non-wetting phase averages
	RecvBuffer.resize(BLOB_AVG_COUNT,NumberComponents_NWP);
/*	for (int b=0; b<NumberComponents_NWP; b++){
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allreduce(&ComponentAverages_NWP(0,b),&RecvBuffer(0),BLOB_AVG_COUNT,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		for (int idx=0; idx<BLOB_AVG_COUNT; idx++) ComponentAverages_NWP(idx,b)=RecvBuffer(idx);
	}
	*/
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&ComponentAverages_NWP(0,0),&RecvBuffer(0,0),BLOB_AVG_COUNT*NumberComponents_NWP,
					MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

	for (int b=0; b<NumberComponents_NWP; b++){
		for (int idx=0; idx<BLOB_AVG_COUNT; idx++) ComponentAverages_NWP(idx,b)=RecvBuffer(b,idx);
	}
	for (int b=0; b<NumberComponents_NWP; b++){
		if (ComponentAverages_NWP(VOL,b) > 0.0){
			double Vn,pn,awn,ans,Jwn,Kwn,lwns,cwns,vsq;

			Vn = ComponentAverages_NWP(VOL,b);
			awn = ComponentAverages_NWP(AWN,b);
			ans = ComponentAverages_NWP(ANS,b);
			van(0) = ComponentAverages_NWP(VX,b)/Vn;
			van(1) = ComponentAverages_NWP(VY,b)/Vn;
			van(2) = ComponentAverages_NWP(VZ,b)/Vn;
			vsq = ComponentAverages_NWP(VSQ,b)/Vn;

			if (ComponentAverages_NWP(TRIMVOL,b) > 0.0){
				pn = ComponentAverages_NWP(PRS,b)/ComponentAverages_NWP(TRIMVOL,b);
			}
			else pn = 0.0;

			if (awn != 0.0){
				Jwn = ComponentAverages_NWP(JWN,b)/awn;
				Kwn = ComponentAverages_NWP(KWN,b)/awn;
				vawn(0) = ComponentAverages_NWP(VWNSX,b)/awn;
				vawn(1) = ComponentAverages_NWP(VWNSY,b)/awn;
				vawn(2) = ComponentAverages_NWP(VWNSZ,b)/awn;
				Gwn(0) = ComponentAverages_NWP(GWNXX,b)/awn;
				Gwn(1) = ComponentAverages_NWP(GWNYY,b)/awn;
				Gwn(2) = ComponentAverages_NWP(GWNZZ,b)/awn;
				Gwn(3) = ComponentAverages_NWP(GWNXY,b)/awn;
				Gwn(4) = ComponentAverages_NWP(GWNXZ,b)/awn;
				Gwn(5) = ComponentAverages_NWP(GWNYZ,b)/awn;
			}
			else Jwn=Kwn=0.0;

			trawn = ComponentAverages_NWP(TRAWN,b);
			trJwn = ComponentAverages_NWP(TRJWN,b);
			if (trawn > 0.0) trJwn /= trawn;

			lwns = ComponentAverages_NWP(LWNS,b);
			if (lwns != 0.0){
				cwns = ComponentAverages_NWP(CWNS,b)/lwns;
				vawns(0) = ComponentAverages_NWP(VWNSX,b)/lwns;
				vawns(1) = ComponentAverages_NWP(VWNSY,b)/lwns;
				vawns(2) = ComponentAverages_NWP(VWNSZ,b)/lwns;
			}
			else  cwns=0.0;

			ComponentAverages_NWP(PRS,b) = pn;
			ComponentAverages_NWP(VX,b) = van(0);
			ComponentAverages_NWP(VY,b) = van(1);
			ComponentAverages_NWP(VZ,b) = van(2);
			ComponentAverages_NWP(VSQ,b) = vsq;

			ComponentAverages_NWP(JWN,b) = Jwn;
			ComponentAverages_NWP(KWN,b) = Kwn;
			ComponentAverages_NWP(VWNX,b) = vawn(0);
			ComponentAverages_NWP(VWNY,b) = vawn(1);
			ComponentAverages_NWP(VWNZ,b) = vawn(2);
			
			ComponentAverages_NWP(GWNXX,b) = Gwn(0);
			ComponentAverages_NWP(GWNYY,b) = Gwn(1);
			ComponentAverages_NWP(GWNZZ,b) = Gwn(2);
			ComponentAverages_NWP(GWNXY,b) = Gwn(3);
			ComponentAverages_NWP(GWNXZ,b) = Gwn(4);
			ComponentAverages_NWP(GWNYZ,b) = Gwn(5);

			ComponentAverages_NWP(CWNS,b) = cwns;
			ComponentAverages_NWP(VWNSX,b) = vawns(0);
			ComponentAverages_NWP(VWNSY,b) = vawns(1);
			ComponentAverages_NWP(VWNSZ,b) = vawns(2);

			ComponentAverages_NWP(CMX,b) /= Vn;
			ComponentAverages_NWP(CMY,b) /= Vn;
			ComponentAverages_NWP(CMZ,b) /= Vn;

			ComponentAverages_NWP(TRJWN,b) = trJwn;

		}
	}

	// reduce the wetting phase averages
	for (int b=0; b<NumberComponents_WP; b++){
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allreduce(&ComponentAverages_WP(0,b),&RecvBuffer(0),BLOB_AVG_COUNT,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		for (int idx=0; idx<BLOB_AVG_COUNT; idx++) ComponentAverages_WP(idx,b)=RecvBuffer(idx);
	}
	
	for (int b=0; b<NumberComponents_WP; b++){
		if (ComponentAverages_WP(VOL,b) > 0.0){
			double Vw,pw,awn,ans,Jwn,Kwn,lwns,cwns,vsq;
			Vw = ComponentAverages_WP(VOL,b);
			awn = ComponentAverages_WP(AWN,b);
			ans = ComponentAverages_WP(ANS,b);
			vaw(0) = ComponentAverages_WP(VX,b)/Vw;
			vaw(1) = ComponentAverages_WP(VY,b)/Vw;
			vaw(2) = ComponentAverages_WP(VZ,b)/Vw;
			vsq = ComponentAverages_WP(VSQ,b)/Vw;

			if (ComponentAverages_WP(TRIMVOL,b) > 0.0){
				pw = ComponentAverages_WP(PRS,b)/ComponentAverages_WP(TRIMVOL,b);
			}
			else pw = 0.0;

			if (awn != 0.0){
				Jwn = ComponentAverages_WP(JWN,b)/awn;
				Kwn = ComponentAverages_WP(KWN,b)/awn;
				vawn(0) = ComponentAverages_WP(VWNSX,b)/awn;
				vawn(1) = ComponentAverages_WP(VWNSY,b)/awn;
				vawn(2) = ComponentAverages_WP(VWNSZ,b)/awn;
				Gwn(0) = ComponentAverages_WP(GWNXX,b)/awn;
				Gwn(1) = ComponentAverages_WP(GWNYY,b)/awn;
				Gwn(2) = ComponentAverages_WP(GWNZZ,b)/awn;
				Gwn(3) = ComponentAverages_WP(GWNXY,b)/awn;
				Gwn(4) = ComponentAverages_WP(GWNXZ,b)/awn;
				Gwn(5) = ComponentAverages_WP(GWNYZ,b)/awn;
			}
			else Jwn=Kwn=0.0;

			trawn = ComponentAverages_WP(TRAWN,b);
			trJwn = ComponentAverages_WP(TRJWN,b);
			if (trawn > 0.0) trJwn /= trawn;

			lwns = ComponentAverages_WP(LWNS,b);
			if (lwns != 0.0){
				cwns = ComponentAverages_WP(CWNS,b)/lwns;
				vawns(0) = ComponentAverages_WP(VWNSX,b)/lwns;
				vawns(1) = ComponentAverages_WP(VWNSY,b)/lwns;
				vawns(2) = ComponentAverages_WP(VWNSZ,b)/lwns;
			}
			else  cwns=0.0;

			ComponentAverages_WP(PRS,b) = pw;
			ComponentAverages_WP(VX,b) = vaw(0);
			ComponentAverages_WP(VY,b) = vaw(1);
			ComponentAverages_WP(VZ,b) = vaw(2);
			ComponentAverages_WP(VSQ,b) = vsq;

			ComponentAverages_WP(JWN,b) = Jwn;
			ComponentAverages_WP(KWN,b) = Kwn;
			ComponentAverages_WP(VWNX,b) = vawn(0);
			ComponentAverages_WP(VWNY,b) = vawn(1);
			ComponentAverages_WP(VWNZ,b) = vawn(2);
			
			ComponentAverages_WP(GWNXX,b) = Gwn(0);
			ComponentAverages_WP(GWNYY,b) = Gwn(1);
			ComponentAverages_WP(GWNZZ,b) = Gwn(2);
			ComponentAverages_WP(GWNXY,b) = Gwn(3);
			ComponentAverages_WP(GWNXZ,b) = Gwn(4);
			ComponentAverages_WP(GWNYZ,b) = Gwn(5);

			ComponentAverages_WP(CWNS,b) = cwns;
			ComponentAverages_WP(VWNSX,b) = vawns(0);
			ComponentAverages_WP(VWNSY,b) = vawns(1);
			ComponentAverages_WP(VWNSZ,b) = vawns(2);

			ComponentAverages_WP(TRJWN,b) = trJwn;
		}
	}
}

void TwoPhase::WriteSurfaces(int logcount){

	int i,j,k,n;
	int ncubes=(Nx-1)*(Ny-1)*(Nz-1);
	Point P,A,B,C;

	std::shared_ptr<IO::TriList> wn_mesh( new IO::TriList() );
	wn_mesh->A.reserve(8*ncubes);
	wn_mesh->B.reserve(8*ncubes);
	wn_mesh->C.reserve(8*ncubes);

	std::shared_ptr<IO::TriList> ns_mesh( new IO::TriList() );
	ns_mesh->A.reserve(8*ncubes);
	ns_mesh->B.reserve(8*ncubes);
	ns_mesh->C.reserve(8*ncubes);

	std::shared_ptr<IO::TriList> ws_mesh( new IO::TriList() );
	ws_mesh->A.reserve(8*ncubes);
	ws_mesh->B.reserve(8*ncubes);
	ws_mesh->C.reserve(8*ncubes);

	for (k=1; k<Nz-1; k++){
		for (j=1; j<Ny-1; j++){
			for (i=1; i<Nx-1; i++){		// Get cube from the list

				//...........................................................................
				// Construct the interfaces and common curve
				pmmc_ConstructLocalCube(SDs, SDn, solid_isovalue, fluid_isovalue,
						nw_pts, nw_tris, Values, ns_pts, ns_tris, ws_pts, ws_tris,
						local_nws_pts, nws_pts, nws_seg, local_sol_pts, local_sol_tris,
						n_local_sol_tris, n_local_sol_pts, n_nw_pts, n_nw_tris,
						n_ws_pts, n_ws_tris, n_ns_tris, n_ns_pts, n_local_nws_pts, n_nws_pts, n_nws_seg,
						i, j, k, Nx, Ny, Nz);

				//.......................................................................................
				// Write the triangle lists to text file
				for (int r=0;r<n_nw_tris;r++){
					A = nw_pts(nw_tris(0,r));
					B = nw_pts(nw_tris(1,r));
					C = nw_pts(nw_tris(2,r));
					// compare the trianlge orientation against the color gradient
					// Orientation of the triangle
					double tri_normal_x = (A.y-B.y)*(B.z-C.z) - (A.z-B.z)*(B.y-C.y);
					double tri_normal_y = (A.z-B.z)*(B.x-C.x) - (A.x-B.x)*(B.z-C.z);
					double tri_normal_z = (A.x-B.x)*(B.y-C.y) - (A.y-B.y)*(B.x-C.x);

					double normal_x = SDn_x(i,j,k);
					double normal_y = SDn_y(i,j,k);
					double normal_z = SDn_z(i,j,k);

					// If the normals don't point in the same direction, flip the orientation of the triangle
					// Right hand rule for triangle orientation is used to determine rendering for most software
					if (normal_x*tri_normal_x + normal_y*tri_normal_y + normal_z*tri_normal_z < 0.0){
						P = A;
						A = C;
						C = P;
					}
					// Remap the points
					A.x += 1.0*Dm.iproc*(Nx-2);
					A.y += 1.0*Dm.jproc*(Nx-2);
					A.z += 1.0*Dm.kproc*(Nx-2);
					B.x += 1.0*Dm.iproc*(Nx-2);
					B.y += 1.0*Dm.jproc*(Nx-2);
					B.z += 1.0*Dm.kproc*(Nx-2);
					C.x += 1.0*Dm.iproc*(Nx-2);
					C.y += 1.0*Dm.jproc*(Nx-2);
					C.z += 1.0*Dm.kproc*(Nx-2);
					wn_mesh->A.push_back(A);
					wn_mesh->B.push_back(B);
					wn_mesh->C.push_back(C);
				}
				for (int r=0;r<n_ws_tris;r++){
					A = ws_pts(ws_tris(0,r));
					B = ws_pts(ws_tris(1,r));
					C = ws_pts(ws_tris(2,r));
					// Remap the points
					A.x += 1.0*Dm.iproc*(Nx-2);
					A.y += 1.0*Dm.jproc*(Nx-2);
					A.z += 1.0*Dm.kproc*(Nx-2);
					B.x += 1.0*Dm.iproc*(Nx-2);
					B.y += 1.0*Dm.jproc*(Nx-2);
					B.z += 1.0*Dm.kproc*(Nx-2);
					C.x += 1.0*Dm.iproc*(Nx-2);
					C.y += 1.0*Dm.jproc*(Nx-2);
					C.z += 1.0*Dm.kproc*(Nx-2);
					ws_mesh->A.push_back(A);
					ws_mesh->B.push_back(B);
					ws_mesh->C.push_back(C);
				}
				for (int r=0;r<n_ns_tris;r++){
					A = ns_pts(ns_tris(0,r));
					B = ns_pts(ns_tris(1,r));
					C = ns_pts(ns_tris(2,r));
					// Remap the points
					A.x += 1.0*Dm.iproc*(Nx-2);
					A.y += 1.0*Dm.jproc*(Nx-2);
					A.z += 1.0*Dm.kproc*(Nx-2);
					B.x += 1.0*Dm.iproc*(Nx-2);
					B.y += 1.0*Dm.jproc*(Nx-2);
					B.z += 1.0*Dm.kproc*(Nx-2);
					C.x += 1.0*Dm.iproc*(Nx-2);
					C.y += 1.0*Dm.jproc*(Nx-2);
					C.z += 1.0*Dm.kproc*(Nx-2);
					ns_mesh->A.push_back(A);
					ns_mesh->B.push_back(B);
					ns_mesh->C.push_back(C);
				}
			}
		}
	}

	std::vector<IO::MeshDataStruct> meshData(4);
	meshData[0].meshName = "wn-tris";
	meshData[0].mesh = wn_mesh;
	meshData[1].meshName = "ws-tris";
	meshData[1].mesh = ws_mesh;
	meshData[2].meshName = "ns-tris";
	meshData[2].mesh = ns_mesh;
	IO::writeData( logcount, meshData, 2);

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
	}
	if (wp_volume_global > 0.0){
		vaw_global(0) = vaw_global(0) / wp_volume_global;
		vaw_global(1) = vaw_global(1) / wp_volume_global;
		vaw_global(2) = vaw_global(2) / wp_volume_global;

	}
	if (vol_n_global > 0.0){
		pan_global = pan_global / vol_n_global;
	}
	if (nwp_volume_global > 0.0){
		van_global(0) = van_global(0) / nwp_volume_global;
		van_global(1) = van_global(1) / nwp_volume_global;
		van_global(2) = van_global(2) / nwp_volume_global;
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

void TwoPhase::PrintComponents(int timestep){
	if (Dm.rank==0){
		printf("PRINT COMPONENT AVEREAGES: time = %i \n",timestep);
		for (int b=0; b<NumberComponents_NWP; b++){
			if (ComponentAverages_NWP(TRIMVOL,b) > 0.0){
				fprintf(NWPLOG,"%i ",timestep-5);
				fprintf(NWPLOG,"%i ",b);
				fprintf(NWPLOG,"%.5g ",ComponentAverages_NWP(VOL,b));
				//			fprintf(NWPLOG,"%.5g ",ComponentAverages_NWP(TRIMVOL,b));
				fprintf(NWPLOG,"%.5g ",ComponentAverages_NWP(PRS,b));
				fprintf(NWPLOG,"%.5g ",ComponentAverages_NWP(AWN,b));
				fprintf(NWPLOG,"%.5g ",ComponentAverages_NWP(ANS,b));
				fprintf(NWPLOG,"%.5g ",ComponentAverages_NWP(JWN,b));
				fprintf(NWPLOG,"%.5g ",ComponentAverages_NWP(KWN,b));
				fprintf(NWPLOG,"%.5g ",ComponentAverages_NWP(LWNS,b));
				fprintf(NWPLOG,"%.5g ",ComponentAverages_NWP(CWNS,b));
				fprintf(NWPLOG,"%.5g ",ComponentAverages_NWP(VX,b));
				fprintf(NWPLOG,"%.5g ",ComponentAverages_NWP(VY,b));
				fprintf(NWPLOG,"%.5g ",ComponentAverages_NWP(VZ,b));
				fprintf(NWPLOG,"%.5g ",ComponentAverages_NWP(VWNX,b));
				fprintf(NWPLOG,"%.5g ",ComponentAverages_NWP(VWNY,b));
				fprintf(NWPLOG,"%.5g ",ComponentAverages_NWP(VWNZ,b));
				fprintf(NWPLOG,"%.5g ",ComponentAverages_NWP(VWNSX,b));
				fprintf(NWPLOG,"%.5g ",ComponentAverages_NWP(VWNSY,b));
				fprintf(NWPLOG,"%.5g ",ComponentAverages_NWP(VWNSZ,b));
				fprintf(NWPLOG,"%.5g ",ComponentAverages_NWP(VSQ,b));
				fprintf(NWPLOG,"%.5g ",ComponentAverages_NWP(GWNXX,b));
				fprintf(NWPLOG,"%.5g ",ComponentAverages_NWP(GWNYY,b));
				fprintf(NWPLOG,"%.5g ",ComponentAverages_NWP(GWNZZ,b));
				fprintf(NWPLOG,"%.5g ",ComponentAverages_NWP(GWNXY,b));
				fprintf(NWPLOG,"%.5g ",ComponentAverages_NWP(GWNXZ,b));
				fprintf(NWPLOG,"%.5g ",ComponentAverages_WP(GWNYZ,b));
				fprintf(NWPLOG,"%.5g ",ComponentAverages_WP(TRAWN,b));
				fprintf(NWPLOG,"%.5g\n",ComponentAverages_WP(TRJWN,b));			}
		}
		fflush(NWPLOG);

		for (int b=0; b<NumberComponents_WP; b++){
			if (ComponentAverages_WP(TRIMVOL,b) > 0.0){
				fprintf(WPLOG,"%i ",timestep-5);
				fprintf(WPLOG,"%i ",b);
				fprintf(WPLOG,"%.5g ",ComponentAverages_WP(VOL,b));
				//			fprintf(WPLOG,"%.5g ",ComponentAverages_WP(TRIMVOL,b));
				fprintf(WPLOG,"%.5g ",ComponentAverages_WP(PRS,b));
				fprintf(WPLOG,"%.5g ",ComponentAverages_WP(AWN,b));
				fprintf(WPLOG,"%.5g ",ComponentAverages_WP(AWS,b));
				fprintf(WPLOG,"%.5g ",ComponentAverages_WP(JWN,b));
				fprintf(WPLOG,"%.5g ",ComponentAverages_WP(KWN,b));
				fprintf(WPLOG,"%.5g ",ComponentAverages_WP(LWNS,b));
				fprintf(WPLOG,"%.5g ",ComponentAverages_WP(CWNS,b));
				fprintf(WPLOG,"%.5g ",ComponentAverages_WP(VX,b));
				fprintf(WPLOG,"%.5g ",ComponentAverages_WP(VY,b));
				fprintf(WPLOG,"%.5g ",ComponentAverages_WP(VZ,b));
				fprintf(WPLOG,"%.5g ",ComponentAverages_WP(VWNX,b));
				fprintf(WPLOG,"%.5g ",ComponentAverages_WP(VWNY,b));
				fprintf(WPLOG,"%.5g ",ComponentAverages_WP(VWNZ,b));
				fprintf(WPLOG,"%.5g ",ComponentAverages_WP(VWNSX,b));
				fprintf(WPLOG,"%.5g ",ComponentAverages_WP(VWNSY,b));
				fprintf(WPLOG,"%.5g ",ComponentAverages_WP(VWNSZ,b));
				fprintf(WPLOG,"%.5g ",ComponentAverages_WP(VSQ,b));
				fprintf(WPLOG,"%.5g ",ComponentAverages_WP(GWNXX,b));
				fprintf(WPLOG,"%.5g ",ComponentAverages_WP(GWNYY,b));
				fprintf(WPLOG,"%.5g ",ComponentAverages_WP(GWNZZ,b));
				fprintf(WPLOG,"%.5g ",ComponentAverages_WP(GWNXY,b));
				fprintf(WPLOG,"%.5g ",ComponentAverages_WP(GWNXZ,b));
				fprintf(WPLOG,"%.5g ",ComponentAverages_WP(GWNYZ,b));
				fprintf(WPLOG,"%.5g ",ComponentAverages_WP(TRAWN,b));
				fprintf(WPLOG,"%.5g\n",ComponentAverages_WP(TRJWN,b));

			}
		}
		fflush(WPLOG);

	}
}

inline int TwoPhase::GetCubeLabel(int i, int j, int k, IntArray &BlobLabel){
	int label;
	label=BlobLabel(i,j,k);
	label=max(label,BlobLabel(i+1,j,k));
	label=max(label,BlobLabel(i,j+1,k));
	label=max(label,BlobLabel(i+1,j+1,k));
	label=max(label,BlobLabel(i,j,k+1));
	label=max(label,BlobLabel(i+1,j,k+1));
	label=max(label,BlobLabel(i,j+1,k+1));
	label=max(label,BlobLabel(i+1,j+1,k+1));
		
	return label;
}

void TwoPhase::SortBlobs(){
  //printf("Sorting the blobs based on volume \n");
  //printf("-----------------------------------------------\n");
	int TempLabel,a,aa,bb,i,j,k,idx;
	double TempValue;
	//.......................................................................
	// Sort NWP components by volume
	//.......................................................................
	IntArray OldLabel(NumberComponents_NWP);
	for (a=0; a<NumberComponents_NWP; a++)	OldLabel(a) = a;
	// Sort the blob averages based on volume
	for (aa=0; aa<NumberComponents_NWP-1; aa++){
		for ( bb=aa+1; bb<NumberComponents_NWP; bb++){
			if (ComponentAverages_NWP(VOL,aa) < ComponentAverages_NWP(VOL,bb)){
				// Exchange location of blobs aa and bb
				//printf("Switch blob %i with %i \n", OldLabel(aa),OldLabel(bb));
				// switch the label
				TempLabel = OldLabel(bb);
				OldLabel(bb) = OldLabel(aa);
				OldLabel(aa) = TempLabel;
				// switch the averages
				for (idx=0; idx<BLOB_AVG_COUNT; idx++){
					TempValue = ComponentAverages_NWP(idx,bb);
					ComponentAverages_NWP(idx,bb) = ComponentAverages_NWP(idx,aa);
					ComponentAverages_NWP(idx,aa) = TempValue;
				}
			}
		}		
	}
	IntArray NewLabel(NumberComponents_NWP);
	for (aa=0; aa<NumberComponents_NWP; aa++){
		// Match the new label for original blob aa
		bb=0;
		while (OldLabel(bb) != aa)	bb++;
		NewLabel(aa) = bb;
	}
	// Re-label the blob ID
	//	printf("Re-labeling the blobs, now indexed by volume \n");
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				if (Label_NWP(i,j,k) > -1){
					TempLabel = NewLabel(Label_NWP(i,j,k));
					Label_NWP(i,j,k) = TempLabel;
				}
			}
		}
	}
	//.......................................................................
	// Sort WP components by volume
	//.......................................................................
	OldLabel.resize(NumberComponents_WP);
	for (a=0; a<NumberComponents_WP; a++)	OldLabel(a) = a;
	// Sort the blob averages based on volume
	for (aa=0; aa<NumberComponents_WP-1; aa++){
		for ( bb=aa+1; bb<NumberComponents_WP; bb++){
			if (ComponentAverages_WP(VOL,aa) < ComponentAverages_WP(VOL,bb)){
				// Exchange location of blobs aa and bb
				//printf("Switch blob %i with %i \n", OldLabel(aa),OldLabel(bb));
				// switch the label
				TempLabel = OldLabel(bb);
				OldLabel(bb) = OldLabel(aa);
				OldLabel(aa) = TempLabel;
				// switch the averages
				for (idx=0; idx<BLOB_AVG_COUNT; idx++){
					TempValue = ComponentAverages_WP(idx,bb);
					ComponentAverages_WP(idx,bb) = ComponentAverages_WP(idx,aa);
					ComponentAverages_WP(idx,aa) = TempValue;
				}
			}
		}		
	}
	NewLabel.resize(NumberComponents_WP);
	for (aa=0; aa<NumberComponents_WP; aa++){
		// Match the new label for original blob aa
		bb=0;
		while (OldLabel(bb) != aa)	bb++;
		NewLabel(aa) = bb;
	}
	// Re-label the blob ID
	//	printf("Re-labeling the blobs, now indexed by volume \n");
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				if (Label_WP(i,j,k) > -1){
					TempLabel = NewLabel(Label_WP(i,j,k));
					Label_WP(i,j,k) = TempLabel;
				}
			}
		}
	}
	//.......................................................................

}

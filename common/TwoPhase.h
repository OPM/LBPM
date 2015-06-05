// Header file for two-phase averaging class
#include "pmmc.h"
#include "Domain.h"
#include "Communication.h"
#include "analysis/analysis.h"
#include <vector>

#define BLOB_AVG_COUNT 26

struct BlobContainer{

    BlobContainer(){
      NBLOBS=0;
    }
    ~BlobContainer(){
    }
    void Set(int size){
      if (NBLOBS!=0) delete [] Data;
      NBLOBS=size;
      Data = new double [BLOB_AVG_COUNT*size];
      //Data.resize(size*BLOB_AVG_COUNT);
      //for (int i=0; i<size*BLOB_AVG_COUNT; i++) Data[i] = 0.0;
    }
    int NBLOBS;
    double * Data;
  //std::vector<double> Data;
    
    // if modified -- make sure to adjust COUNT so that
    // there is enough memory to save all the averages
    double Vn(int IDX){return Data[BLOB_AVG_COUNT*IDX];}
    double pan(int IDX){return Data[BLOB_AVG_COUNT*IDX+1];}
    double awn(int IDX){return Data[BLOB_AVG_COUNT*IDX+2];}
    double ans(int IDX){return Data[BLOB_AVG_COUNT*IDX+3];}
    double Jwn(int IDX){return Data[BLOB_AVG_COUNT*IDX+4];}
    double Kwn(int IDX){return Data[BLOB_AVG_COUNT*IDX+5];}
    double lwns(int IDX){return Data[BLOB_AVG_COUNT*IDX+6];}
    double cwns(int IDX){return Data[BLOB_AVG_COUNT*IDX+7];}
    double vanx(int IDX){return Data[BLOB_AVG_COUNT*IDX+8];}
    double vany(int IDX){return Data[BLOB_AVG_COUNT*IDX+9];}
    double vanz(int IDX){return Data[BLOB_AVG_COUNT*IDX+10];}
    double vawnx(int IDX){return Data[BLOB_AVG_COUNT*IDX+11];}
    double vawny(int IDX){return Data[BLOB_AVG_COUNT*IDX+12];}
    double vawnz(int IDX){return Data[BLOB_AVG_COUNT*IDX+13];}
    double Gwnxx(int IDX){return Data[BLOB_AVG_COUNT*IDX+14];}
    double Gwnyy(int IDX){return Data[BLOB_AVG_COUNT*IDX+15];}
    double Gwnzz(int IDX){return Data[BLOB_AVG_COUNT*IDX+16];}
    double Gwnxy(int IDX){return Data[BLOB_AVG_COUNT*IDX+17];}
    double Gwnxz(int IDX){return Data[BLOB_AVG_COUNT*IDX+18];}
    double Gwnyz(int IDX){return Data[BLOB_AVG_COUNT*IDX+19];}
    double Gnsxx(int IDX){return Data[BLOB_AVG_COUNT*IDX+20];}
    double Gnsyy(int IDX){return Data[BLOB_AVG_COUNT*IDX+22];}
    double Gnszz(int IDX){return Data[BLOB_AVG_COUNT*IDX+23];}
    double Gnsxy(int IDX){return Data[BLOB_AVG_COUNT*IDX+23];}
    double Gnsxz(int IDX){return Data[BLOB_AVG_COUNT*IDX+24];}
    double Gnsyz(int IDX){return Data[BLOB_AVG_COUNT*IDX+25];}

    void Vn(int IDX, double value){ Data[BLOB_AVG_COUNT*IDX]+=value;}
    void pan(int IDX, double value){ Data[BLOB_AVG_COUNT*IDX+1]+=value;}
    void awn(int IDX, double value){ Data[BLOB_AVG_COUNT*IDX+2]+=value;}
    void ans(int IDX, double value){ Data[BLOB_AVG_COUNT*IDX+3]+=value;}
    void Jwn(int IDX, double value){ Data[BLOB_AVG_COUNT*IDX+4]+=value;}
    void Kwn(int IDX, double value){ Data[BLOB_AVG_COUNT*IDX+5]+=value;}
    void lwns(int IDX, double value){ Data[BLOB_AVG_COUNT*IDX+6]+=value;}
    void cwns(int IDX, double value){ Data[BLOB_AVG_COUNT*IDX+7]+=value;}
    void vanx(int IDX, double value){ Data[BLOB_AVG_COUNT*IDX+8]+=value;}
    void vany(int IDX, double value){ Data[BLOB_AVG_COUNT*IDX+9]+=value;}
    void vanz(int IDX, double value){ Data[BLOB_AVG_COUNT*IDX+10]+=value;}
    void vawnx(int IDX, double value){ Data[BLOB_AVG_COUNT*IDX+11]+=value;}
    void vawny(int IDX, double value){ Data[BLOB_AVG_COUNT*IDX+12]+=value;}
    void vawnz(int IDX, double value){ Data[BLOB_AVG_COUNT*IDX+13]+=value;}
    void Gwnxx(int IDX, double value){ Data[BLOB_AVG_COUNT*IDX+14]+=value;}
    void Gwnyy(int IDX, double value){ Data[BLOB_AVG_COUNT*IDX+15]+=value;}
    void Gwnzz(int IDX, double value){ Data[BLOB_AVG_COUNT*IDX+16]+=value;}
    void Gwnxy(int IDX, double value){ Data[BLOB_AVG_COUNT*IDX+17]+=value;}
    void Gwnxz(int IDX, double value){ Data[BLOB_AVG_COUNT*IDX+18]+=value;}
    void Gwnyz(int IDX, double value){ Data[BLOB_AVG_COUNT*IDX+19]+=value;}
    void Gnsxx(int IDX, double value){ Data[BLOB_AVG_COUNT*IDX+20]+=value;}
    void Gnsyy(int IDX, double value){ Data[BLOB_AVG_COUNT*IDX+22]+=value;}
    void Gnszz(int IDX, double value){ Data[BLOB_AVG_COUNT*IDX+23]+=value;}
    void Gnsxy(int IDX, double value){ Data[BLOB_AVG_COUNT*IDX+23]+=value;}
    void Gnsxz(int IDX, double value){ Data[BLOB_AVG_COUNT*IDX+24]+=value;}
    void Gnsyz(int IDX, double value){ Data[BLOB_AVG_COUNT*IDX+25]+=value;}

};

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

public:
	Domain& Dm;
	int nblobs_global;
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
	IntArray BlobLabel;
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
	//	BlobContainer BlobAverages;
	DoubleArray BlobAverages;
	//...........................................................................
	TwoPhase(Domain &dm) : Dm(dm){
		Nx=dm.Nx; Ny=dm.Ny; Nz=dm.Nz;
		Volume=(Nx-2)*(Ny-2)*(Nz-2)*Dm.nprocx*Dm.nprocy*Dm.nprocz*1.0;

		ncubes=(Nx-2)*(Ny-2)*(Nz-2);
		cubeList.resize(3,ncubes);

		// Global arrays
		BlobLabel.resize(Nx,Ny,Nz);
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
	}
	~TwoPhase(){

	}

	void Initialize();
	void SetupCubes(Domain &Dm);
	void UpdateMeshValues();
	void UpdateSolid();
	void ComputeDelPhi();
	void ColorToSignedDistance(double Beta, double *ColorData, double *DistData);
	void ComputeLocal();
	void ComputeLocalBlob();
	void Reduce();
	void NonDimensionalize(double D, double viscosity, double IFT);
	void PrintAll(int timestep);
	int GetCubeLabel(int i, int j, int k);
	void SortBlobs();

};

void TwoPhase::ColorToSignedDistance(double Beta, double *ColorData, double *DistData){

	double factor,temp,value;
	factor=0.5/Beta;
/*	for (int n=0; n<Nx*Ny*Nz; n++){
	 	value = ColorData[n];
		if (value > 0.999 ) DistData[n] = 4.0;
		else if (value < -0.999 ) DistData[n] = -4.0;
		else 	DistData[n] = factor*log((1.0+value)/(1.0-value));
	}

	// Initialize to -1,1 (segmentation)
	for (int k=0; k<Nz; k++){
		for (int j=0; j<Ny; j++){
			for (int i=0; i<Nx; i++){
				value = ColorData[n];
				temp = factor*log((1.0+value)/(1.0-value));
				if (temp > 1.0) DistData(i,j,k) = 1.0;
				else if (temp < -1.0) DistData(i,j,k) = -1.0;
				else DistData(i,j,k) = temp;
			}
		}
	}
*/
	for (int n=0; n<Nx*Ny*Nz; n++)	DistData[n] = ColorData[n];
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

}
void TwoPhase::ComputeLocal(){
	int i,j,k,n;
	int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};

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
					if (DelPhi(n) < 1e-4){
						vol_n += 0.125;
						// pressure
						pan += 0.125*Press(n);
						// velocity
						van(0) += 0.125*Vel_x(n);
						van(1) += 0.125*Vel_y(n);
						van(2) += 0.125*Vel_z(n);
					}
				}
				else{
					wp_volume += 0.125;
					if (DelPhi(n) < 1e-4){
						// volume the excludes the interfacial region
						vol_w += 0.125;
						// pressure
						paw += 0.125*Press(n);
						// velocity
						vaw(0) += 0.125*Vel_x(n);
						vaw(1) += 0.125*Vel_y(n);
						vaw(2) += 0.125*Vel_z(n);
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

			pmmc_CurveCurvature(SDn, SDs, KNwns_values, KGwns_values, KNwns, KGwns,
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

void TwoPhase::ComputeLocalBlob(){
    int i,j,k,n,label;
	double vF,vS;
	vF = 0.0; vS=0.0;
//    const RankInfoStruct rank_info(Dm.rank,Dm.nprocx,Dm.nprocy,Dm.nprocz);

	int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};
        // get the maximum label locally -- then compute number of global blobs
	label=0; 
	nblobs_global = 0;
	for (n=0; n<Nx*Ny*Nz; n++){
	  if (label < BlobLabel(n)) label = BlobLabel(n);
	}
	MPI_Allreduce(&label,&nblobs_global,1,MPI_INT,MPI_MAX,Dm.Comm);
	nblobs_global+=1;
	if (Dm.rank==0) printf("Number of blobs is %i \n",nblobs_global);

	//BlobAverages.Set(nblobs_global);
	BlobAverages.resize(BLOB_AVG_COUNT,nblobs_global);
    BlobAverages.fill(0.0);
	// Perform averaging
	for (int c=0;c<ncubes;c++){
		// Get cube from the list
		i = cubeList(0,c);
		j = cubeList(1,c);
		k = cubeList(2,c);
		label=GetCubeLabel(i,j,k);

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
					BlobAverages(1,label) += 0.125;
					// volume the excludes the interfacial region
					if (DelPhi(n) < 1e-4){
						BlobAverages(0,label) += 0.125;
						// pressure
						BlobAverages(2,label ) += 0.125*Press(n);
						// velocity
						BlobAverages(9,label) += 0.125*Vel_x(n);
						BlobAverages(10,label) += 0.125*Vel_y(n);
						BlobAverages(11,label) += 0.125*Vel_z(n);
					}
				}

				else{
					wp_volume += 0.125;
					if (DelPhi(n) < 1e-4){
						// volume the excludes the interfacial region
						vol_w += 0.125;
						// pressure
						if (isnan(Press(n))) printf("Pressure is nan!\n");
						else paw += 0.125*Press(n);
						// velocity
						vaw(0) += 0.125*Vel_x(n);
						vaw(1) += 0.125*Vel_y(n);
						vaw(2) += 0.125*Vel_z(n);
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
		if (n_nw_pts> 0){
			BlobAverages(9,label) += pmmc_CubeSurfaceInterpValue(CubeValues,MeanCurvature,nw_pts,nw_tris,Values,i,j,k,n_nw_pts,n_nw_tris);
			BlobAverages(5,label) += pmmc_CubeSurfaceInterpValue(CubeValues,GaussCurvature,nw_pts,nw_tris,Values,i,j,k,n_nw_pts,n_nw_tris);

			// Integrate the trimmed mean curvature (hard-coded to use a distance of 4 pixels)
			pmmc_CubeTrimSurfaceInterpValues(CubeValues,MeanCurvature,SDs,nw_pts,nw_tris,Values,DistanceValues,
					i,j,k,n_nw_pts,n_nw_tris,trimdist,BlobAverages(12,label),BlobAverages(13,label));

			pmmc_CubeTrimSurfaceInterpInverseValues(CubeValues,MeanCurvature,SDs,nw_pts,nw_tris,Values,DistanceValues,
					i,j,k,n_nw_pts,n_nw_tris,trimdist,dummy,trRwn);

			// Compute the normal speed of the interface
			pmmc_InterfaceSpeed(dPdt, SDn_x, SDn_y, SDn_z, CubeValues, nw_pts, nw_tris,
					NormalVector, InterfaceSpeed, vawn, i, j, k, n_nw_pts, n_nw_tris);

			BlobAverages(3,label) += pmmc_CubeSurfaceOrientation(Gwn,nw_pts,nw_tris,n_nw_tris);

		}

		// Common curve averages
		if (n_local_nws_pts > 0){
			BlobAverages(8,label) += pmmc_CubeContactAngle(CubeValues,Values,SDn_x,SDn_y,SDn_z,SDs_x,SDs_y,SDs_z,
					local_nws_pts,i,j,k,n_local_nws_pts);

			pmmc_CommonCurveSpeed(CubeValues, dPdt, vawns, SDn_x, SDn_y, SDn_z,SDs_x,SDs_y,SDs_z,
					local_nws_pts,i,j,k,n_local_nws_pts);

			pmmc_CurveCurvature(SDn, SDs, KNwns_values, KGwns_values, KNwns, KGwns,
					nws_pts, n_nws_pts, i, j, k);

			BlobAverages(7,label) +=  pmmc_CubeCurveLength(local_nws_pts,n_local_nws_pts);
		}

		// Solid interface averages
		if (n_local_sol_pts > 0){
			As  += pmmc_CubeSurfaceArea(local_sol_pts,local_sol_tris,n_local_sol_tris);
			// Compute the surface orientation and the interfacial area
			BlobAverages(4,label) += pmmc_CubeSurfaceOrientation(Gns,ns_pts,ns_tris,n_ns_tris);
			aws += pmmc_CubeSurfaceOrientation(Gws,ws_pts,ws_tris,n_ws_tris);
			//...........................................................................
		}
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

inline int TwoPhase::GetCubeLabel(int i, int j, int k){
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
	IntArray OldLabel(nblobs_global);
	for (a=0; a<nblobs_global; a++)	OldLabel(a) = a;
	// Sort the blob averages based on volume
	for (aa=0; aa<nblobs_global-1; aa++){
		for ( bb=aa+1; bb<nblobs_global; bb++){
			if (BlobAverages(0,aa) < BlobAverages(0,bb)){
				// Exchange location of blobs aa and bb
				//printf("Switch blob %i with %i \n", OldLabel(aa),OldLabel(bb));
				// switch the label
				TempLabel = OldLabel(bb);
				OldLabel(bb) = OldLabel(aa);
				OldLabel(aa) = TempLabel;
				// switch the averages
				for (idx=0; idx<BLOB_AVG_COUNT; idx++){
					TempValue = BlobAverages(idx,bb);
					BlobAverages(idx,bb) = BlobAverages(idx,aa);
					BlobAverages(idx,aa) = TempValue;
				}
			}
		}		
	}
	
	IntArray NewLabel(nblobs_global);
	for (aa=0; aa<nblobs_global; aa++){
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
				if (BlobLabel(i,j,k) > -1){
					TempLabel = NewLabel(BlobLabel(i,j,k));
					BlobLabel(i,j,k) = TempLabel;
				}
			}
		}
	}
}

// Header file for two-phase averaging class
#include "pmmc.h"
#include "Domain.h"

class TwoPhase{

	//...........................................................................
	int n_nw_pts,n_ns_pts,n_ws_pts,n_nws_pts,n_local_sol_pts,n_local_nws_pts;
	int n_nw_tris,n_ns_tris,n_ws_tris,n_nws_seg,n_local_sol_tris;
	int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};
	//...........................................................................
	int nc;
	int kstart,kfinish;

	double fluid_isovalue, solid_isovalue;
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

public:
	int ncubes;
	//...........................................................................
	// Averaging variables
	//...........................................................................
	// local averages (to each MPI process)
	double trimdist; 						// pixel distance to trim surface for specified averages
	double awn,ans,aws,lwns,nwp_volume;
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
	double nwp_volume_global;					// volume for the wetting phase (for saturation)
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
	DoubleArray Phase_x;		// Gradient of the phase indicator field
	DoubleArray Phase_y;
	DoubleArray Phase_z;
	DoubleArray Phase_tplus;
	DoubleArray Phase_tminus;
	DoubleArray Vel_x;		// Gradient of the phase indicator field
	DoubleArray Vel_y;
	DoubleArray Vel_z;
	//...........................................................................
	TwoPhase(Domain &dm){
		Nx=dm.Nx; Ny=dm.Ny; Nz=dm.Nz;

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
		Phase_x.New(Nx,Ny,Nz);		// Gradient of the phase indicator field
		Phase_y.New(Nx,Ny,Nz);
		Phase_z.New(Nx,Ny,Nz);
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

	}
	~TwoPhase(){
	}

	void Initialize();
	void SetupCubes(Domain &Dm);
	void UpdateMeshValues();
	void PhaseToSignedDistance(double Beta);
	void ComputeLocal();
	void Reduce();
	void PrintAll(int timestep, FILE*);

};

void TwoPhase::PhaseToSignedDistance(double Beta){

	double temp=0.5/Beta;
	for (int n=0; n<Nx*Ny*Nz; n++){
	  double value = Phase.data[n];
		SDn.data[n] = temp*log((1.0+value)/(1.0-value));
	}
}

void TwoPhase::Initialize(){

	trimdist=1.0;
	fluid_isovalue=solid_isovalue=0.0;
	// Initialize the averaged quantities
	awn = aws = ans = lwns = 0.0;
	nwp_volume = 0.0;
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
	if (Dm.pBC && Dm.kproc==0)				kstart = 4;
	if (Dm.pBC && Dm.kproc==Dm.nprocz-1)	kfinish = Nz-4;

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

void TwoPhase::UpdateMeshValues(){

	//...........................................................................
	// Compute the gradients of the phase indicator and signed distance fields
	pmmc_MeshGradient(Phase,Phase_x,Phase_y,Phase_z,Nx,Ny,Nz);
	pmmc_MeshGradient(SDs,SDs_x,SDs_y,SDs_z,Nx,Ny,Nz);
	//...........................................................................
	// Compute the mesh curvature of the phase indicator field
	pmmc_MeshCurvature(SDn, MeanCurvature, GaussCurvature, Nx, Ny, Nz);
	//...........................................................................
	// Update the time derivative of non-dimensional density field
	for (int n=0; n<Nx*Ny*Nz; n++)	dPdt(n) = 0.1*(Phase_tplus(n) - Phase_tminus(n));
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
				delphi = sqrt(Phase_x.data[n]*Phase_x.data[n]+Phase_y.data[n]*Phase_y.data[n]+Phase_z.data[n]*Phase_z.data[n]);
				// Compute the non-wetting phase volume contribution
				if ( Phase(i+cube[p][0],j+cube[p][1],k+cube[p][2]) > 0 ){
					nwp_volume += 0.125;
					// volume the excludes the interfacial region
					if (delphi < 1e-4){
						vol_n += 0.125;
						// pressure
						pan += 0.125*Press.data[n];
						// velocity
						van(0) += 0.125*Vel_x.data[n];
						van(1) += 0.125*Vel_y.data[n];
						van(2) += 0.125*Vel_z.data[n];
					}
				}
				else if (delphi < 1e-4){
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

		//...........................................................................
		// Construct the interfaces and common curve
		pmmc_ConstructLocalCube(SDs, SDn, solid_isovalue, fluid_isovalue,
				nw_pts, nw_tris, Values, ns_pts, ns_tris, ws_pts, ws_tris,
				local_nws_pts, nws_pts, nws_seg, local_sol_pts, local_sol_tris,
				n_local_sol_tris, n_local_sol_pts, n_nw_pts, n_nw_tris,
				n_ws_pts, n_ws_tris, n_ns_tris, n_ns_pts, n_local_nws_pts, n_nws_pts, n_nws_seg,
				i, j, k, Nx, Ny, Nz);

		// Integrate the contact angle
		efawns += pmmc_CubeContactAngle(CubeValues,Values,Phase_x,Phase_y,Phase_z,SDs_x,SDs_y,SDs_z,
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
		pmmc_InterfaceSpeed(dPdt, Phase_x, Phase_y, Phase_z, CubeValues, nw_pts, nw_tris,
				NormalVector, InterfaceSpeed, vawn, i, j, k, n_nw_pts, n_nw_tris);

		pmmc_CommonCurveSpeed(CubeValues, dPdt, vawns,Phase_x,Phase_y,Phase_z,SDs_x,SDs_y,SDs_z,
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
		//printf("finished cube %i \n",c);
	}
}

void TwoPhase::Reduce(){
  awn_global=awn;
  Jwn_global=Jwn/awn;
}

void TwoPhase::PrintAll(int timestep, FILE *TIMELOG){
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


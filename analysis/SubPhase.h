/*
 * Sub-phase averaging tools
 */

class SubPhase{
public:
	double Volume;
	// input variables
	double rho_n, rho_w;
	double nu_n, nu_w;
	double gamma_wn;
	double Fx, Fy, Fz;
	
	// mass
	double Mwc,Mwd,Mwi,Mnc,Mnd,Mni;
	// momentum
	double Pwc_x,Pwd_x,Pwi_x,Pnc_x,Pnd_x,Pni_x;
	double Pwc_y,Pwd_y,Pwi_y,Pnc_y,Pnd_y,Pni_y;
	double Pwc_z,Pwd_z,Pwi_z,Pnc_z,Pnd_z,Pni_z;
	// energy
	double Kwc,Kwd,Kwi,Knc,Knd,Kni;

	// Geometric measures
	double Vwc, Awc, Hwc, Xwc; 	// connected w
	double Vwd, Awd, Hwd, Xwd; 	// disconnected w
	double Vnc, Anc, Hnc, Xnc; 	// connected n
	double Vnd, And, Hnd, Xnd; 	// disconnected n
	double Vi, Ai, Hi, Xi;		// interface

	// other measures
	double pnc,pnd,pni,pwc,pwd,pwi;
	
	//...........................................................................
    int Nx,Ny,Nz;
	IntArray PhaseID;		// Phase ID array (solid=0, non-wetting=1, wetting=2)
	BlobIDArray Label_WP;   // Wetting phase label
	BlobIDArray Label_NWP;  // Non-wetting phase label index (0:nblobs-1)
	std::vector<BlobIDType> Label_NWP_map;  // Non-wetting phase label for each index
	DoubleArray Density;	// density field 
	DoubleArray Phi;		// phase indicator field
	DoubleArray DelPhi;		// Magnitude of Gradient of the phase indicator field
	DoubleArray Pressure; 	// pressure field
	DoubleArray Vel_x;		// velocity field
	DoubleArray Vel_y;
	DoubleArray Vel_z;

	std::shared_ptr<Minkowski> morph_w;
	std::shared_ptr<Minkowski> morph_n;
	std::shared_ptr<Minkowski> morph_i;

	SubPhase(std::shared_ptr <Domain> Dm);
	~SubPhase();
	
	void SetParams(double rhoA, double rhoB, double tauA, double tauB, double force_x, double force_y, double force_z, double alpha);
	PhaseAverages();
	
private:
};

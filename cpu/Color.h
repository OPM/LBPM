extern "C" void dvc_InitDenColor(char *ID, double *Den, double *Phi, double das, double dbs, int Nx, int Ny, int Nz, int S);
extern "C" void dvc_InitDenColorDistance(char *ID, double *Den, double *Phi, double *Distance,
								double das, double dbs, double beta, double xp, int Nx, int Ny, int Nz, int S);
extern "C" void dvc_Compute_VELOCITY(char *ID, double *disteven, double *distodd, double *vel, int Nx, int Ny, int Nz, int S);
extern "C" void dvc_ComputePressureD3Q19(char *ID, double *disteven, double *distodd, double *Pressure, 
									int Nx, int Ny, int Nz, int S);
extern "C" void dvc_PressureBC_inlet(double *disteven, double *distodd, double din,
								  int Nx, int Ny, int Nz, int S);
extern "C" void dvc_PressureBC_outlet(double *disteven, double *distodd, double dout,
								   int Nx, int Ny, int Nz, int S, int outlet);
extern "C" void dvc_ComputeColorGradient(char *ID, double *phi, double *ColorGrad, int Nx, int Ny, int Nz, int S);
extern "C" void dvc_ColorCollide( char *ID, double *disteven, double *distodd, double *ColorGrad,
								double *Velocity, int Nx, int Ny, int Nz, int S,double rlx_setA, double rlx_setB,
								double alpha, double beta, double Fx, double Fy, double Fz, bool pBC);
extern "C" void dvc_ColorCollideOpt( char *ID, double *disteven, double *distodd, double *phi, double *ColorGrad,
								double *Velocity, int Nx, int Ny, int Nz, int S,double rlx_setA, double rlx_setB, 
								double alpha, double beta, double Fx, double Fy, double Fz);
extern "C" void dvc_DensityStreamD3Q7(char *ID, double *Den, double *Copy, double *Phi, double *ColorGrad, double *Velocity,
		double beta, int Nx, int Ny, int Nz, bool pBC, int S);
extern "C" void dvc_ComputePhi(char *ID, double *Phi, double *Copy, double *Den, int N, int S);

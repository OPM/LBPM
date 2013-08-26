extern void InitDenColor(char *ID, double *Den, double *Phi, double das, double dbs, int N);
extern void InitD3Q19(char *ID, double *f_even, double *f_odd, int Nx, int Ny, int Nz);

extern void Compute_VELOCITY(char *ID, double *disteven, double *distodd, double *vel, int Nx, int Ny, int Nz);

//*************************************************************************
//*************************************************************************
extern  void PressureBC_inlet(double *disteven, double *distodd, double din,
			      int Nx, int Ny, int Nz);
extern  void PressureBC_outlet(double *disteven, double *distodd, double dout,
			       int Nx, int Ny, int Nz, int S, int outlet);
//*************************************************************************
extern void ComputeColorGradient(char *ID, double *phi, double *ColorGrad, int Nx, int Ny, int Nz);
//*************************************************************************
extern void ColorCollide( char *ID, double *disteven, double *distodd, double *ColorGrad,
						 double *Velocity, int Nx, int Ny, int Nz, double rlx_setA, double rlx_setB,
			  double alpha, double beta, double Fx, double Fy, double Fz, bool pBC);
//*************************************************************************
extern void DensityStreamD3Q7(char *ID, double *Den, double *Copy, double *Phi, double *ColorGrad, double *Velocity,
			      double beta, int Nx, int Ny, int Nz, bool pBC);
extern void ComputePhi(char *ID, double *Phi, double *Copy, double *Den, int N);

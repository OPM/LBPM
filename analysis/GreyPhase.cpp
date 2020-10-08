#include "analysis/SubPhase.h"

// Constructor
GreyPhase::GreyPhase(std::shared_ptr <Domain> dm):
	Dm(dm)
{
	Nx=dm->Nx; Ny=dm->Ny; Nz=dm->Nz;
	Volume=(Nx-2)*(Ny-2)*(Nz-2)*Dm->nprocx()*Dm->nprocy()*Dm->nprocz()*1.0;
	
	// Global arrays
	PhaseID.resize(Nx,Ny,Nz);       PhaseID.fill(0);
	Rho_n.resize(Nx,Ny,Nz);       	Rho_n.fill(0);
	Rho_w.resize(Nx,Ny,Nz);       	Rho_w.fill(0);
	Pressure.resize(Nx,Ny,Nz);      Pressure.fill(0);
	Phi.resize(Nx,Ny,Nz);         	Phi.fill(0);
	DelPhi.resize(Nx,Ny,Nz);        DelPhi.fill(0);
	Vel_x.resize(Nx,Ny,Nz);         Vel_x.fill(0);	    // Gradient of the phase indicator field
	Vel_y.resize(Nx,Ny,Nz);         Vel_y.fill(0);
	Vel_z.resize(Nx,Ny,Nz);         Vel_z.fill(0);
	//.........................................
	
	if (Dm->rank()==0){
		bool WriteHeader=false;
		TIMELOG = fopen("timelog.csv","r");
		if (TIMELOG != NULL)
			fclose(TIMELOG);
		else
			WriteHeader=true;

		TIMELOG = fopen("timelog.csv","a+");
		if (WriteHeader)
		{
			// If timelog is empty, write a short header to list the averages
			//fprintf(TIMELOG,"--------------------------------------------------------------------------------------\n");
			fprintf(TIMELOG,"sw krw krn vw vn pw pn\n");				
		}
	}
}


// Destructor
GreyPhase::~GreyPhase()
{

}

void GreyPhase::Write(int timestep)
{

}

void GreyPhase::SetParams(double rhoA, double rhoB, double tauA, double tauB, double force_x, double force_y, double force_z, double alpha, double B)
{
	Fx = force_x;
	Fy = force_y;
	Fz = force_z;
	rho_n = rhoA;
	rho_w = rhoB;
	nu_n = (tauA-0.5)/3.f;
	nu_w = (tauB-0.5)/3.f;
	gamma_wn = 5.796*alpha;
	beta = B;
}

void GreyPhase::Basic(){
	int i,j,k,n,imin,jmin,kmin,kmax;

	// If external boundary conditions are set, do not average over the inlet
	kmin=1; kmax=Nz-1;
	imin=jmin=1;

	double count_w = 0.0;
	double count_n = 0.0;
	
	for (k=0; k<Nz; k++){
		for (j=0; j<Ny; j++){
			for (i=0; i<Nx; i++){
				n = k*Nx*Ny + j*Nx + i;
				// Compute volume averages
				if ( Dm->id[n] > 0 ){
					// compute density
					double nA = Rho_n(n);
					double nB = Rho_w(n);
					double phi = (nA-nB)/(nA+nB);
					Phi(n) = phi;
				}
			}
		}
	}
	
	for (k=kmin; k<kmax; k++){
		for (j=jmin; j<Ny-1; j++){
			for (i=imin; i<Nx-1; i++){
				n = k*Nx*Ny + j*Nx + i;
				// Compute volume averages
				if ( Dm->id[n] > 0 ){
					// compute density
					double nA = Rho_n(n);
					double nB = Rho_w(n);
					double phi = (nA-nB)/(nA+nB);					
					/*					if ( phi > 0.0 ){
						nA = 1.0;
						nb.V += 1.0;
						nb.M += nA*rho_n;						
						// velocity
						nb.Px += rho_n*nA*Vel_x(n);
						nb.Py += rho_n*nA*Vel_y(n);
						nb.Pz += rho_n*nA*Vel_z(n);
					}
					else{
						nB = 1.0;
						wb.M += nB*rho_w;
						wb.V += 1.0;

						// velocity
						wb.Px += rho_w*nB*Vel_x(n);
						wb.Py += rho_w*nB*Vel_y(n);
						wb.Pz += rho_w*nB*Vel_z(n);
					}
					if ( phi > 0.99 ){
						nb.p += Pressure(n);
						count_n += 1.0;
					}
					else if ( phi < -0.99 ){
						wb.p += Pressure(n);
						count_w += 1.0;
					}
					*/
				}
			}
		}
	}
	/*	gwb.V=sumReduce( Dm->Comm, wb.V);
	gnb.V=sumReduce( Dm->Comm, nb.V);
	gwb.M=sumReduce( Dm->Comm, wb.M);
	gnb.M=sumReduce( Dm->Comm, nb.M);
	gwb.Px=sumReduce( Dm->Comm, wb.Px);
	gwb.Py=sumReduce( Dm->Comm, wb.Py);
	gwb.Pz=sumReduce( Dm->Comm, wb.Pz);
	gnb.Px=sumReduce( Dm->Comm, nb.Px);
	gnb.Py=sumReduce( Dm->Comm, nb.Py);
	gnb.Pz=sumReduce( Dm->Comm, nb.Pz);
	
	count_w=sumReduce( Dm->Comm, count_w);
	count_n=sumReduce( Dm->Comm, count_n);
	if (count_w > 0.0)
		gwb.p=sumReduce( Dm->Comm, wb.p) / count_w;
	else 
		gwb.p = 0.0;
	if (count_n > 0.0)
		gnb.p=sumReduce( Dm->Comm, nb.p) / count_n;
	else 
		gnb.p = 0.0;

	// check for NaN
	bool err=false;
	if (gwb.V != gwb.V) err=true;
	if (gnb.V != gnb.V) err=true;
	if (gwb.p != gwb.p) err=true;
	if (gnb.p != gnb.p) err=true;
	if (gwb.Px != gwb.Px) err=true;
	if (gwb.Py != gwb.Py) err=true;
	if (gwb.Pz != gwb.Pz) err=true;
	if (gnb.Px != gnb.Px) err=true;
	if (gnb.Py != gnb.Py) err=true;
	if (gnb.Pz != gnb.Pz) err=true;	
	
	if (Dm->rank() == 0){
		double force_mag = sqrt(Fx*Fx+Fy*Fy+Fz*Fz);
		double dir_x = 0.0;
		double dir_y = 0.0;
		double dir_z = 0.0;
		if (force_mag > 0.0){
			dir_x = Fx/force_mag;
			dir_y = Fy/force_mag;
			dir_z = Fz/force_mag;
		}
		else {
			// default to z direction
			dir_x = 0.0;
			dir_y = 0.0;
			dir_z = 1.0;
		}
		if (Dm->BoundaryCondition == 1 || Dm->BoundaryCondition == 2 || Dm->BoundaryCondition == 3 || Dm->BoundaryCondition == 4 ){
			// compute the pressure drop
			double pressure_drop = (Pressure(Nx*Ny + Nx + 1) - 1.0) / 3.0;
			double length = ((Nz-2)*Dm->nprocz());
			force_mag -= pressure_drop/length;
		}
		if (force_mag == 0.0){
			// default to z direction
			dir_x = 0.0;
			dir_y = 0.0;
			dir_z = 1.0;
			force_mag = 1.0;
		}
		double saturation=gwb.V/(gwb.V + gnb.V);
		double water_flow_rate=gwb.V*(gwb.Px*dir_x + gwb.Py*dir_y + gwb.Pz*dir_z)/gwb.M / Dm->Volume;
		double not_water_flow_rate=gnb.V*(gnb.Px*dir_x + gnb.Py*dir_y + gnb.Pz*dir_z)/gnb.M/ Dm->Volume;
		//double total_flow_rate = water_flow_rate + not_water_flow_rate;
		//double fractional_flow = water_flow_rate / total_flow_rate;

		double h = Dm->voxel_length;		
		double krn = h*h*nu_n*not_water_flow_rate / force_mag ;
		double krw = h*h*nu_w*water_flow_rate / force_mag;
		//printf("   water saturation = %f, fractional flow =%f \n",saturation,fractional_flow);
		fprintf(TIMELOG,"%.5g %.5g %.5g %.5g %.5g %.5g %.5g\n",saturation,krw,krn,h*water_flow_rate,h*not_water_flow_rate, gwb.p, gnb.p); 
		fflush(TIMELOG); 
	}
*/
	if (err==true){
		// exception if simulation produceds NaN
		printf("GreyPhase.cpp: NaN encountered, may need to check simulation parameters \n");
	}	
	ASSERT(err==false);
}
/*
inline void InterfaceTransportMeasures( double beta, double rA, double rB, double nA, double nB, 
		double nx, double ny, double nz, double ux, double uy, double uz, interface &I){
	
	double A1,A2,A3,A4,A5,A6;
	double B1,B2,B3,B4,B5,B6;
	double nAB,delta;
	// Instantiate mass transport distributions
	// Stationary value - distribution 0
	nAB = 1.0/(nA+nB);
	//...............................................
	// q = 0,2,4
	// Cq = {1,0,0}, {0,1,0}, {0,0,1}
	delta = beta*nA*nB*nAB*0.1111111111111111*nx;
	if (!(nA*nB*nAB>0)) delta=0;
	A1 = nA*(0.1111111111111111*(1+4.5*ux))+delta;
	B1 = nB*(0.1111111111111111*(1+4.5*ux))-delta;
	A2 = nA*(0.1111111111111111*(1-4.5*ux))-delta;
	B2 = nB*(0.1111111111111111*(1-4.5*ux))+delta;

	//...............................................
	// Cq = {0,1,0}
	delta = beta*nA*nB*nAB*0.1111111111111111*ny;
	if (!(nA*nB*nAB>0)) delta=0;
	A3 = nA*(0.1111111111111111*(1+4.5*uy))+delta;
	B3 = nB*(0.1111111111111111*(1+4.5*uy))-delta;
	A4 = nA*(0.1111111111111111*(1-4.5*uy))-delta;
	B4 = nB*(0.1111111111111111*(1-4.5*uy))+delta;

	//...............................................
	// q = 4
	// Cq = {0,0,1}
	delta = beta*nA*nB*nAB*0.1111111111111111*nz;
	if (!(nA*nB*nAB>0)) delta=0;
	A5 = nA*(0.1111111111111111*(1+4.5*uz))+delta;
	B5 = nB*(0.1111111111111111*(1+4.5*uz))-delta;
	A6 = nA*(0.1111111111111111*(1-4.5*uz))-delta;
	B6 = nB*(0.1111111111111111*(1-4.5*uz))+delta;
	
	double unx = (A1-A2);
	double uny = (A3-A4);
	double unz = (A5-A6);
	double uwx = (B1-B2);
	double uwy = (B3-B4);
	double uwz = (B5-B6);
	
	I.Mn += rA*nA;
	I.Mw += rB*nB;
	I.Pnx += rA*nA*unx;
	I.Pny += rA*nA*uny;
	I.Pnz += rA*nA*unz;
	I.Pwx += rB*nB*uwx;
	I.Pwy += rB*nB*uwy;
	I.Pwz += rB*nB*uwz;
	I.Kn += rA*nA*(unx*unx + uny*uny + unz*unz);
	I.Kw += rB*nB*(uwx*uwx + uwy*uwy + uwz*uwz);

}
*/

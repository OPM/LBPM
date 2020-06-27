#include "analysis/distance.h"



/******************************************************************
* A fast distance calculation                                     *
******************************************************************/
template<class TYPE>
void CalcDist( Array<TYPE> &Distance, const Array<char> &ID, const Domain &Dm,
    const std::array<bool,3>& periodic, const std::array<double,3>& dx )
{
    ASSERT( Distance.size() == ID.size() );
    std::array<int,3> n = { Dm.Nx-2, Dm.Ny-2, Dm.Nz-2 };
    fillHalo<int> fillData(  Dm.Comm, Dm.rank_info, n, {1,1,1}, 50, 1, {true,false,false}, periodic );
    Array<int> id(ID.size());
    Array<Vec> vecDist(Distance.size());
    for (size_t i=0; i<ID.length(); i++)
        id(i) = ID(i) == 0 ? -1:1;
    fillData.fill( id );
    CalcVecDist( vecDist, id, Dm, periodic, dx );
    for (size_t i=0; i<Distance.length(); i++)
        Distance(i) = id(i)*vecDist(i).norm();
}


/******************************************************************
* Vector-based distance calculation                               *
* Initialize cells adjacent to boundaries                         *
******************************************************************/
static void calcVecInitialize( Array<Vec> &d, const Array<int> &ID, double dx, double dy, double dz )
{
    d.fill( Vec( 1e50, 1e50, 1e50 ) );
    const double dx0 = 0.5*dx;
    const double dy0 = 0.5*dy;
    const double dz0 = 0.5*dz;
    //const double dxy0 = 0.25*sqrt( dx*dx + dy*dy );
    //const double dxz0 = 0.25*sqrt( dx*dx + dz*dz );
    //const double dyz0 = 0.25*sqrt( dy*dy + dz*dz );
    //const double dxyz0 = sqrt( dx*dx + dy*dy + dz*dz );
    int Nx = d.size(0);
    int Ny = d.size(1);
    int Nz = d.size(2);
    for (int k=1; k<Nz-1; k++) {
        for (int j=1; j<Ny-1; j++) {
            for (int i=1; i<Nx-1; i++) {
                int id = ID(i,j,k);
                bool x[2] = { id != ID(i-1,j,k), id != ID(i+1,j,k) };
                bool y[2] = { id != ID(i,j-1,k), id != ID(i,j+1,k) };
                bool z[2] = { id != ID(i,j,k-1), id != ID(i,j,k+1) };
                if ( x[0] )  d(i,j,k) = Vec( dx0, 0, 0 );
                if ( x[1] )  d(i,j,k) = Vec( -dx0, 0, 0 );
                if ( y[0] )  d(i,j,k) = Vec( 0, dy0, 0 );
                if ( y[1] )  d(i,j,k) = Vec( 0, -dy0, 0 );
                if ( z[0] )  d(i,j,k) = Vec( 0, 0, dz0 );
                if ( z[1] )  d(i,j,k) = Vec( 0, 0, -dz0 );
                /*if ( x[0] && y[0] )  d(i,j,k) = Vec( dxy0, dxy0, 0 );
                if ( x[0] && y[1] )  d(i,j,k) = Vec( dxy0, -dxy0, 0 );
                if ( x[1] && y[0] )  d(i,j,k) = Vec( -dxy0, dxy0, 0 );
                if ( x[1] && y[1] )  d(i,j,k) = Vec( -dxy0, -dxy0, 0 );
                if ( x[0] && z[0] )  d(i,j,k) = Vec( dxz0, 0, dxz0 ); 
                if ( x[0] && z[1] )  d(i,j,k) = Vec( dxz0, 0, -dxz0 );
                if ( x[1] && z[0] )  d(i,j,k) = Vec( -dxz0, 0, dxz0 );
                if ( x[1] && z[1] )  d(i,j,k) = Vec( -dxz0, 0, -dxz0 );
                if ( y[0] && z[0] )  d(i,j,k) = Vec( 0, dyz0, dyz0 );
                if ( y[0] && z[1] )  d(i,j,k) = Vec( 0, dyz0, -dyz0 );
                if ( y[1] && z[0] )  d(i,j,k) = Vec( 0, -dyz0, dyz0 );
                if ( y[1] && z[1] )  d(i,j,k) = Vec( 0, -dyz0, -dyz0 );*/
            }
        }
    }

}


/******************************************************************
* Vector-based distance calculation                               *
* Update interior cells                                           *
******************************************************************/
static double calcVecUpdateInterior( Array<Vec> &d, double dx, double dy, double dz )
{
    double err = 0;
    int Nx = d.size(0);
    int Ny = d.size(1);
    int Nz = d.size(2);
    // Propagate (+,+,+)
    for (int k=1; k<Nz; k++) {
        for (int j=1; j<Ny; j++) {
            for (int i=1; i<Nx; i++) {
                auto vx = d(i-1,j,k);
                auto vy = d(i,j-1,k);
                auto vz = d(i,j,k-1);
                vx.x += dx;
                vy.y += dy;
                vz.z += dz;
                auto v = std::min( std::min(vx,vy), vz );
                double d1 = v.norm2();
                double d2 = d(i,j,k).norm2();
                if ( d1 < d2 ) {
                    d(i,j,k) = v;
                    err = std::max( err, sqrt(d2)-sqrt(d1) );
                }
            }
        }
    }
    // Propagate (-,-,-)
    for (int k=Nz-2; k>=0; k--) {
        for (int j=Ny-2; j>=0; j--) {
            for (int i=Nx-2; i>=0; i--) {
                auto vx = d(i+1,j,k);
                auto vy = d(i,j+1,k);
                auto vz = d(i,j,k+1);
                vx.x -= dx;
                vy.y -= dy;
                vz.z -= dz;
                auto v = std::min( std::min(vx,vy), vz );
                double d1 = v.norm2();
                double d2 = d(i,j,k).norm2();
                if ( d1 < d2 ) {
                    d(i,j,k) = v;
                    err = std::max( err, sqrt(d2)-sqrt(d1) );
                }
            }
        }
    }
    return err;
}


/******************************************************************
* Vector-based distance calculation                               *
******************************************************************/
void CalcVecDist( Array<Vec> &d, const Array<int> &ID0, const Domain &Dm,
    const std::array<bool,3>& periodic, const std::array<double,3>& dx )
{
    std::array<int,3> N = { Dm.Nx, Dm.Ny, Dm.Nz };
    std::array<int,3> n = { Dm.Nx-2, Dm.Ny-2, Dm.Nz-2 };
    // Create ID with ghosts
    Array<int> ID(N[0],N[1],N[2]);
    fillHalo<int> fillDataID(  Dm.Comm, Dm.rank_info, n, {1,1,1}, 50, 1, {true,true,true}, periodic );
    fillDataID.copy( ID0, ID );
    // Fill ghosts with nearest neighbor
    for (int k=1; k<N[2]-1; k++) {
        for (int j=1; j<N[1]-1; j++) {
            ID(0,j,k) = ID(1,j,k);
            ID(N[0]-1,j,k) = ID(N[0]-2,j,k);
        }
    }
    for (int k=1; k<N[2]-1; k++) {
        for (int i=0; i<N[0]; i++) {
            ID(i,0,k) = ID(i,1,k);
            ID(i,N[1]-1,k) = ID(i,N[1]-2,k);
        }
    }
    for (int i=0; i<N[0]; i++) {
        for (int j=0; j<N[1]; j++) {
            ID(i,j,0) = ID(i,j,1);
            ID(i,j,N[2]-1) = ID(i,j,N[2]-2);
        }
    }
    // Communicate ghosts
    fillDataID.fill( ID );
    // Create communicator for distance
    fillHalo<Vec> fillData( Dm.Comm, Dm.rank_info, n, {1,1,1}, 50, 1, {true,false,false}, periodic );
    // Calculate the local distances
    calcVecInitialize( d, ID, dx[0], dx[1], dx[2] );
    double err = 1e100;
    double tol = 0.5 * std::min( std::min(dx[0],dx[1]), dx[2] );
    for (int it=0; it<=50 && err>tol; it++) {
        err = calcVecUpdateInterior( d, dx[0], dx[1], dx[2] );
    }
    // Calculate the global distances
    int N_it = Dm.nprocx() + Dm.nprocy() + Dm.nprocz() + 100;
    for ( int it=0; it<N_it; it++ ) {
        // Update ghosts
        fillData.fill( d );
        // Update distance
        double err = calcVecUpdateInterior( d, dx[0], dx[1], dx[2] );
        // Check if we are finished
        err = maxReduce( Dm.Comm, err );
        if ( err < tol )
            break;
    }
}

double Eikonal(DoubleArray &Distance, char *ID, Domain &Dm, int timesteps){

  /*
   * This routine converts the data in the Distance array to a signed distance
   * by solving the equation df/dt = sign(1-|grad f|), where Distance provides
   * the values of f on the mesh associated with domain Dm
   * It has been tested with segmented data initialized to values [-1,1]
   * and will converge toward the signed distance to the surface bounding the associated phases
   *
   * Reference:
   * Min C (2010) On reinitializing level set functions, Journal of Computational Physics229
   */

  int i,j,k;
  double dt=0.1;
  double Dx,Dy,Dz;
  double Dxp,Dxm,Dyp,Dym,Dzp,Dzm;
  double Dxxp,Dxxm,Dyyp,Dyym,Dzzp,Dzzm;
  double sign,norm;
  double LocalVar,GlobalVar,LocalMax,GlobalMax;

  int xdim,ydim,zdim;
  xdim=Dm.Nx-2;
  ydim=Dm.Ny-2;
  zdim=Dm.Nz-2;
  fillHalo<double> fillData(Dm.Comm, Dm.rank_info,xdim,ydim,zdim,1,1,1,0,1);

  // Arrays to store the second derivatives
  DoubleArray Dxx(Dm.Nx,Dm.Ny,Dm.Nz);
  DoubleArray Dyy(Dm.Nx,Dm.Ny,Dm.Nz);
  DoubleArray Dzz(Dm.Nx,Dm.Ny,Dm.Nz);

  int count = 0;
  while (count < timesteps){

    // Communicate the halo of values
    fillData.fill(Distance);

    // Compute second order derivatives
    for (k=1;k<Dm.Nz-1;k++){
      for (j=1;j<Dm.Ny-1;j++){
	for (i=1;i<Dm.Nx-1;i++){
	  Dxx(i,j,k) = Distance(i+1,j,k) + Distance(i-1,j,k) - 2*Distance(i,j,k);
	  Dyy(i,j,k) = Distance(i,j+1,k) + Distance(i,j-1,k) - 2*Distance(i,j,k);
	  Dzz(i,j,k) = Distance(i,j,k+1) + Distance(i,j,k-1) - 2*Distance(i,j,k);
	}
      }
    }
    fillData.fill(Dxx);
    fillData.fill(Dyy);
    fillData.fill(Dzz);

    LocalMax=LocalVar=0.0;
    // Execute the next timestep
    for (k=1;k<Dm.Nz-1;k++){
      for (j=1;j<Dm.Ny-1;j++){
	for (i=1;i<Dm.Nx-1;i++){

	  int n = k*Dm.Nx*Dm.Ny + j*Dm.Nx + i;

	  sign = -1;
	  if (ID[n] == 1) sign = 1;

	  // local second derivative terms
	  Dxxp = minmod(Dxx(i,j,k),Dxx(i+1,j,k));
	  Dyyp = minmod(Dyy(i,j,k),Dyy(i,j+1,k));
	  Dzzp = minmod(Dzz(i,j,k),Dzz(i,j,k+1));
	  Dxxm = minmod(Dxx(i,j,k),Dxx(i-1,j,k));
	  Dyym = minmod(Dyy(i,j,k),Dyy(i,j-1,k));
	  Dzzm = minmod(Dzz(i,j,k),Dzz(i,j,k-1));

	  /* //............Compute upwind derivatives ...................
                    Dxp = Distance(i+1,j,k) - Distance(i,j,k) + 0.5*Dxxp;
                    Dyp = Distance(i,j+1,k) - Distance(i,j,k) + 0.5*Dyyp;
                    Dzp = Distance(i,j,k+1) - Distance(i,j,k) + 0.5*Dzzp;
                    Dxm = Distance(i,j,k) - Distance(i-1,j,k) + 0.5*Dxxm;
                    Dym = Distance(i,j,k) - Distance(i,j-1,k) + 0.5*Dyym;
                    Dzm = Distance(i,j,k) - Distance(i,j,k-1) + 0.5*Dzzm;
	  */
	  Dxp = Distance(i+1,j,k)- Distance(i,j,k) - 0.5*Dxxp;
	  Dyp = Distance(i,j+1,k)- Distance(i,j,k) - 0.5*Dyyp;
	  Dzp = Distance(i,j,k+1)- Distance(i,j,k) - 0.5*Dzzp;

	  Dxm = Distance(i,j,k) - Distance(i-1,j,k) + 0.5*Dxxm;
	  Dym = Distance(i,j,k) - Distance(i,j-1,k) + 0.5*Dyym;
	  Dzm = Distance(i,j,k) - Distance(i,j,k-1) + 0.5*Dzzm;

	  // Compute upwind derivatives for Godunov Hamiltonian
	  if (sign < 0.0){
	    if (Dxp + Dxm > 0.f)  Dx = Dxp*Dxp;
	    else Dx = Dxm*Dxm;

	    if (Dyp + Dym > 0.f)  Dy = Dyp*Dyp;
	    else Dy = Dym*Dym;

	    if (Dzp + Dzm > 0.f)  Dz = Dzp*Dzp;
	    else Dz = Dzm*Dzm;
	  }
	  else{

	    if (Dxp + Dxm < 0.f)  Dx = Dxp*Dxp;
	    else Dx = Dxm*Dxm;

	    if (Dyp + Dym < 0.f)  Dy = Dyp*Dyp;
	    else Dy = Dym*Dym;

	    if (Dzp + Dzm < 0.f)  Dz = Dzp*Dzp;
	    else Dz = Dzm*Dzm;
	  }

	  //Dx = max(Dxp*Dxp,Dxm*Dxm);
	  //Dy = max(Dyp*Dyp,Dym*Dym);
	  //Dz = max(Dzp*Dzp,Dzm*Dzm);

	  norm=sqrt(Dx + Dy + Dz);
	  if (norm > 1.0) norm=1.0;

	  Distance(i,j,k) += dt*sign*(1.0 - norm);
	  LocalVar +=  dt*sign*(1.0 - norm);

	  if (fabs(dt*sign*(1.0 - norm)) > LocalMax)
	    LocalMax = fabs(dt*sign*(1.0 - norm));
	}
      }
    }

    MPI_Allreduce(&LocalVar,&GlobalVar,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
    MPI_Allreduce(&LocalMax,&GlobalMax,1,MPI_DOUBLE,MPI_MAX,Dm.Comm);
    GlobalVar /= (Dm.Nx-2)*(Dm.Ny-2)*(Dm.Nz-2)*Dm.nprocx*Dm.nprocy*Dm.nprocz;
    count++;

    if (count%50 == 0 && Dm.rank==0 )
      printf("Time=%i, Max variation=%f, Global variation=%f \n",count,GlobalMax,GlobalVar);

    if (fabs(GlobalMax) < 1e-5){
      if (Dm.rank==0) printf("Exiting with max tolerance of 1e-5 \n");
      count=timesteps;
    }
  }
  return GlobalVar;
}

// Explicit instantiations
template void CalcDist<float>( Array<float>&, const Array<char>&, const Domain&, const std::array<bool,3>&, const std::array<double,3>& );
template void CalcDist<double>( Array<double>&, const Array<char>&, const Domain&, const std::array<bool,3>&, const std::array<double,3>& );



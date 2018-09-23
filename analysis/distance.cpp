/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University

  This file is part of the Open Porous Media project (OPM).
  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
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


// Explicit instantiations
template void CalcDist<float>( Array<float>&, const Array<char>&, const Domain&, const std::array<bool,3>&, const std::array<double,3>& );
template void CalcDist<double>( Array<double>&, const Array<char>&, const Domain&, const std::array<bool,3>&, const std::array<double,3>& );



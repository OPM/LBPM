#include "analysis/distance.h"


// Check if we need to recompute distance after updating ghose values
static bool checkUpdate( const Array<Vec> &d, double dx, double dy, double dz )
{
    auto s = d.size();
    bool test[3] = { false, false, false };
    // Check x-direction
    Vec v1, v2;
    for (size_t k=1; k<s[2]-1; k++) {
        for (size_t j=1; j<s[1]-1; j++) {
            v1 = d(1,j,k);
            v2 = d(0,j,k);
            v2.x += dx;
            test[0] = test[0] || v2 < v1;
            v1 = d(s[0]-2,j,k);
            v2 = d(s[0]-1,j,k);
            v2.x -= dx;
            test[0] = test[0] || v2 < v1;
        }
    }
    // Check y-direction
    for (size_t k=1; k<s[2]-1; k++) {
        for (size_t i=1; i<s[0]-1; i++) {
            v1 = d(i,1,k);
            v2 = d(i,0,k);
            v2.y += dy;
            test[1] = test[1] || v2 < v1;
            v1 = d(i,s[1]-2,k);
            v2 = d(i,s[1]-1,k);
            v2.y -= dy;
            test[1] = test[1] || v2 < v1;
        }
    }
    // Check z-direction
    for (size_t j=1; j<s[1]-1; j++) {
        for (size_t i=1; i<s[0]-1; i++) {
            v1 = d(i,j,1);
            v2 = d(i,j,0);
            v2.z += dz;
            test[2] = test[2] || v2 < v1;
            v1 = d(i,j,s[2]-2);
            v2 = d(i,j,s[2]-1);
            v2.z -= dz;
            test[2] = test[2] || v2 < v1;
        }
    }
    return test[0] || test[1] || test[2];
}


/******************************************************************
* A fast distance calculation                                     *
******************************************************************/
template<class TYPE>
void CalcDist( Array<TYPE> &Distance, const Array<char> &ID, const Domain &Dm, const std::array<bool,3>& periodic )
{
    ASSERT( Distance.size() == ID.size() );
    std::array<int,3> n = { Dm.Nx-2, Dm.Ny-2, Dm.Nz-2 };
    fillHalo<int> fillData(  Dm.Comm, Dm.rank_info, n, {1,1,1}, 50, 1, {true,false,false}, periodic );
    Array<int> id(ID.size());
    Array<Vec> vecDist(Distance.size());
    for (size_t i=0; i<ID.length(); i++)
        id(i) = ID(i) == 0 ? -1:1;
    fillData.fill( id );
    CalcVecDist( vecDist, id, Dm, periodic );
    for (size_t i=0; i<Distance.length(); i++)
        Distance(i) = id(i)*vecDist(i).norm();
}
void CalcVecDist( Array<Vec> &d, const Array<int> &ID0, const Domain &Dm, const std::array<bool,3>& periodic )
{
    const double dx = 1.0;
    const double dy = 1.0;
    const double dz = 1.0;
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
            ID(i,j,N[1]-1) = ID(i,j,N[2]-2);
        }
    }
    // Communicate ghosts
    fillDataID.fill( ID );
    // Create communicator for distance
    fillHalo<Vec> fillData( Dm.Comm, Dm.rank_info, n, {1,1,1}, 50, 1, {true,false,false}, periodic );
    // Initialize the vector distance
    d.fill( Vec( 1e50, 1e50, 1e50 ) );
    const double dx0 = 0.5*dx;
    const double dy0 = 0.5*dy;
    const double dz0 = 0.5*dz;
    //const double dxy0 = 0.25*sqrt( dx*dx + dy*dy );
    //const double dxz0 = 0.25*sqrt( dx*dx + dz*dz );
    //const double dyz0 = 0.25*sqrt( dy*dy + dz*dz );
    //const double dxyz0 = sqrt( dx*dx + dy*dy + dz*dz );
    for (int k=1; k<Dm.Nz-1; k++) {
        for (int j=1; j<Dm.Ny-1; j++) {
            for (int i=1; i<Dm.Nx-1; i++) {
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
    int N_it = Dm.nprocx() + Dm.nprocy() + Dm.nprocz() + 3;
    for ( int it=0; it<N_it; it++ ) {
        //if ( Dm.rank() == 0 )
        //    printf("Runing iteration %i\n",it+1);
        // Propagate +/- x-direction
        for (int k=0; k<Dm.Nz; k++) {
            for (int j=0; j<Dm.Ny; j++) {
                for (int i=1; i<Dm.Nx; i++) {
                    auto v1 = d(i,j,k);
                    auto v2 = d(i-1,j,k);
                    v2.x += dx;
                    if ( v2 < v1 )
                        d(i,j,k) = v2;
                }
                for (int i=Dm.Nx-2; i>=0; i--) {
                    auto v1 = d(i,j,k);
                    auto v2 = d(i+1,j,k);
                    v2.x -= dx;
                    if ( v2 < v1 )
                        d(i,j,k) = v2;
                }
            }
        }
        // Propagate +/- y-direction
        for (int k=0; k<Dm.Nz; k++) {
            for (int i=0; i<Dm.Nx; i++) {
                for (int j=1; j<Dm.Ny; j++) {
                    auto v1 = d(i,j,k);
                    auto v2 = d(i,j-1,k);
                    v2.y += dy;
                    if ( v2 < v1 )
                        d(i,j,k) = v2;
                }
                for (int j=Dm.Ny-2; j>=0; j--) {
                    auto v1 = d(i,j,k);
                    auto v2 = d(i,j+1,k);
                    v2.y -= dy;
                    if ( v2 < v1 )
                        d(i,j,k) = v2;
                }
            }
        }
        // Propagate +/- z-direction
        for (int j=0; j<Dm.Ny; j++) {
            for (int i=0; i<Dm.Nx; i++) {
                for (int k=1; k<Dm.Nz; k++) {
                    auto v1 = d(i,j,k);
                    auto v2 = d(i,j,k-1);
                    v2.z += dz;
                    if ( v2 < v1 )
                        d(i,j,k) = v2;
                }
                for (int k=Dm.Nz-2; k>=0; k--) {
                    auto v1 = d(i,j,k);
                    auto v2 = d(i,j,k+1);
                    v2.z -= dz;
                    if ( v2 < v1 )
                        d(i,j,k) = v2;
                }
            }
        }
        fillData.fill( d );
        bool test = checkUpdate( d, dx, dy, dz );
        test = sumReduce( Dm.Comm, test );
        if ( !test )
            break;
    }
}
template void CalcDist<float>( Array<float>&, const Array<char>&, const Domain&, const std::array<bool,3>& );
template void CalcDist<double>( Array<double>&, const Array<char>&, const Domain&, const std::array<bool,3>& );



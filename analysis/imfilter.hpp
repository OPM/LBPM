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
#include "analysis/imfilter.h"
#include "ProfilerApp.h"
#include <math.h>
#include <string.h>


#define IMFILTER_INSIST INSIST
#define IMFILTER_ASSERT ASSERT
#define IMFILTER_ERROR  ERROR


// Function to convert an index
static inline int imfilter_index( int index, const int N, const imfilter::BC bc )
{
    if ( index < 0 || index >= N ) {
        if ( bc == imfilter::BC::symmetric ) {
            index = ( 2 * N - index ) % N;
        } else if ( bc == imfilter::BC::replicate ) {
            index = index < 0 ? 0 : N - 1;
        } else if ( bc == imfilter::BC::circular ) {
            index = ( index + N ) % N;
        } else if ( bc == imfilter::BC::fixed ) {
            index = -1;
        }
    }
    return index;
}


// Function to copy a 1D array and pad with the appropriate BC
template<class TYPE>
static inline void copy_array( const int N, const int Ns, const int Nh,
    const TYPE *A, const imfilter::BC BC, const TYPE X, TYPE *B )
{
    // Fill the center with a memcpy
    for (int i=0; i<N; i++ )
        B[i+Nh] = A[i*Ns];
    // Fill the boundaries
    for (int i=0; i<Nh; i++ ) {
        int j1 = imfilter_index( -(i+1), N, BC );
        int j2 = imfilter_index(   N+i, N, BC );
        B[Nh-i-1] = j1==-1 ? X : B[Nh+j1];
        B[N+Nh+i] = j2==-1 ? X : B[Nh+j2];
    }
}


/********************************************************
* Perform a 1D filter in a single direction             *
********************************************************/
template<class TYPE>
static void filter_direction( int Ns, int N, int Ne, int Nh, const TYPE *H,
    imfilter::BC boundary, TYPE X, TYPE *A )
{
    if ( Nh < 0 )
        IMFILTER_ERROR("Invalid filter size");
    if ( Nh == 0 ) {
        for (int i=0; i<Ns*N*Ne; i++)
            A[i] *= H[0];
        return;
    }
    TYPE *tmp = new TYPE[N+2*Nh];
    for (int j=0; j<Ne; j++) {
        for (int i=0; i<Ns; i++) {
            copy_array( N, Ns, Nh, &A[i+j*Ns*N], boundary, X, tmp );
            for (int k=0; k<N; k++) {
                TYPE tmp2 = 0;
                for (int m=0; m<=2*Nh; m++)
                    tmp2 += H[m] * tmp[k+m];
                A[i+k*Ns+j*Ns*N] = tmp2;
            }
        }
    }
    delete[] tmp;
}
template<class TYPE>
static void filter_direction( int Ns, int N, int Ne, int Nh,
    std::function<TYPE(const Array<TYPE>&)> H, imfilter::BC boundary, TYPE X, TYPE *A )
{
    if ( Nh < 0 )
        IMFILTER_ERROR("Invalid filter size");
    TYPE *tmp = new TYPE[N+2*Nh];
    Array<TYPE> tmp2(2*Nh+1);
    for (int j=0; j<Ne; j++) {
        for (int i=0; i<Ns; i++) {
            copy_array( N, Ns, Nh, &A[i+j*Ns*N], boundary, X, tmp );
            for (int k=0; k<N; k++) {
                for (int m=0; m<=2*Nh; m++)
                    tmp2(m) = tmp[k+m];
                A[i+k*Ns+j*Ns*N] = H(tmp2);
            }
        }
    }
    delete[] tmp;
}
template<class TYPE>
static void filter_direction( int Ns, int N, int Ne, int Nh,
    std::function<TYPE(int, const TYPE*)> H, imfilter::BC boundary, TYPE X, TYPE *A )
{
    if ( Nh < 0 )
        IMFILTER_ERROR("Invalid filter size");
    TYPE *tmp = new TYPE[N+2*Nh];
    int Nh2 = 2*Nh+1;
    for (int j=0; j<Ne; j++) {
        for (int i=0; i<Ns; i++) {
            copy_array( N, Ns, Nh, &A[i+j*Ns*N], boundary, X, tmp );
            for (int k=0; k<N; k++)
                A[i+k*Ns+j*Ns*N] = H(Nh2,&tmp[k]);
        }
    }
    delete[] tmp;
}



/********************************************************
*  Create a filter                                      *
********************************************************/
template<class TYPE>
Array<TYPE> imfilter::create_filter( const std::vector<int>& N0, const std::string &type, const void *args )
{
    std::vector<size_t> N2(N0.size());
    for (size_t i=0; i<N2.size(); i++)
        N2[i] = 2*N0[i]+1;
    Array<TYPE> h(N2);
    h.fill(0);
    if ( type == "average" ) {
        // average
        h.fill( 1.0 / static_cast<TYPE>( h.length() ) );
    } else if ( type == "gaussian" ) {
        // gaussian
        if ( N0.size() > 3 )
            IMFILTER_ERROR( "Not implimented for dimensions > 3" );
        TYPE std[3] = { 0.5, 0.5, 0.5 };
        if ( args != NULL ) {
            const TYPE *args2 = reinterpret_cast<const TYPE*>( args );
            for ( size_t d = 0; d < N0.size(); d++ )
                std[d]  = args2[d];
        }
        auto N = N0;
        N.resize(3,0);
        for ( int k = -N[2]; k <= N[2]; k++ ) {
            for ( int j = -N[1]; j <= N[1]; j++ ) {
                for ( int i = -N[0]; i <= N[0]; i++ ) {
                    h(i+N[0],j+N[1],k+N[2]) = 
                        exp( -i * i / ( 2 * std[0] * std[0] ) ) *
                        exp( -j * j / ( 2 * std[1] * std[1] ) ) *
                        exp( -k * k / ( 2 * std[2] * std[2] ) );
                }
            }
        }
        h.scale( 1.0/h.sum() );
    } else {
        IMFILTER_ERROR( "Unknown filter" );
    }
    return h;
}


// Perform 2-D filtering
template<class TYPE>
void imfilter_2D( int Nx, int Ny, const TYPE *A, int Nhx, int Nhy, const TYPE *H,
    imfilter::BC BCx, imfilter::BC BCy, const TYPE X, TYPE *B )
{
    IMFILTER_ASSERT( A != B );
    PROFILE_START( "imfilter_2D" );
    memset( B, 0, Nx * Ny * sizeof( TYPE ) );
    for ( int j1 = 0; j1 < Ny; j1++ ) {
        for ( int i1 = 0; i1 < Nx; i1++ ) {
            TYPE tmp = 0;
            if ( i1 >= Nhx && i1 < Nx - Nhx && j1 >= Nhy && j1 < Ny - Nhy ) {
                int ijkh = 0;
                for ( int j2 = j1 - Nhy; j2 <= j1 + Nhy; j2++ ) {
                    for ( int i2 = i1 - Nhx; i2 <= i1 + Nhx; i2++, ijkh++ )
                        tmp += H[ijkh] * A[i2 + j2 * Nx];
                }
            } else {
                int ijkh = 0;
                for ( int jh = -Nhy; jh <= Nhy; jh++ ) {
                    int j2 = imfilter_index( j1+jh, Ny, BCy );
                    for ( int ih = -Nhx; ih <= Nhx; ih++ ) {
                        int i2     = imfilter_index( i1+ih, Nx, BCx );
                        bool fixed = i2 == -1 || j2 == -1;
                        TYPE A2  = fixed ? X : A[i2 + j2 * Nx];
                        tmp += H[ijkh] * A2;
                        ijkh++;
                    }
                }
            }
            B[i1 + j1 * Nx] = tmp;
        }
    }
    PROFILE_STOP( "imfilter_2D" );
}


// Perform 3-D filtering
template<class TYPE>
void imfilter_3D( int Nx, int Ny, int Nz, const TYPE *A, int Nhx, int Nhy, int Nhz,
    const TYPE *H, imfilter::BC BCx, imfilter::BC BCy, imfilter::BC BCz,
    const TYPE X, TYPE *B )
{
    IMFILTER_ASSERT( A != B );
    PROFILE_START( "imfilter_3D" );
    memset( B, 0, Nx * Ny * Nz * sizeof( TYPE ) );
    for ( int k1 = 0; k1 < Nz; k1++ ) {
        for ( int j1 = 0; j1 < Ny; j1++ ) {
            for ( int i1 = 0; i1 < Nx; i1++ ) {
                TYPE tmp = 0;
                int ijkh   = 0;
                for ( int kh = -Nhz; kh <= Nhz; kh++ ) {
                    int k2 = imfilter_index( k1+kh, Nz, BCz );
                    for ( int jh = -Nhy; jh <= Nhy; jh++ ) {
                        int j2 = imfilter_index( j1+jh, Ny, BCy );
                        for ( int ih = -Nhx; ih <= Nhx; ih++ ) {
                            int i2     = imfilter_index( i1+ih, Nx, BCx );
                            bool fixed = i2 == -1 || j2 == -1 || k2 == -1;
                            TYPE A2  = fixed ? X : A[i2 + j2 * Nx + k2 * Nx * Ny];
                            tmp += H[ijkh] * A2;
                            ijkh++;
                        }
                    }
                }
                B[i1 + j1 * Nx + k1 * Nx * Ny] = tmp;
            }
        }
    }
    PROFILE_STOP( "imfilter_3D" );
}


/********************************************************
*  Perform N-D filtering                                *
********************************************************/
template<class TYPE>
Array<TYPE> imfilter::imfilter( const Array<TYPE>& A, 
    const Array<TYPE>& H, const std::vector<imfilter::BC>& BC, const TYPE X )
{
    IMFILTER_ASSERT( A.ndim() == H.ndim() );
    IMFILTER_ASSERT( A.ndim() == BC.size() );
    std::vector<size_t> Nh = H.size();
    for (int d=0; d<A.ndim(); d++) {
        Nh[d] = (H.size(d)-1)/2;
        IMFILTER_INSIST(2*Nh[d]+1==H.size(d),"Filter must be of size 2*N+1");
    }
    auto B = A;
    if ( A.ndim() == 1 ) {
        PROFILE_START( "imfilter_1D" );
        filter_direction( 1, A.size(0), 1, Nh[0], H.data(), BC[0], X, B.data() );
        PROFILE_STOP( "imfilter_1D" );
    } else if ( A.ndim() == 2 ) {
        imfilter_2D( A.size(0), A.size(1), A.data(), Nh[0], Nh[1], H.data(), BC[0], BC[1], X, B.data() );
    } else if ( A.ndim() == 3 ) {
        imfilter_3D( A.size(0), A.size(1), A.size(2), A.data(),
            Nh[0], Nh[1], Nh[2], H.data(), BC[0], BC[1], BC[2], X, B.data() );
    } else {
        IMFILTER_ERROR( "Arbitrary dimension not yet supported" );
    }
    return B;
}
template<class TYPE>
Array<TYPE> imfilter::imfilter( const Array<TYPE>& A, const std::vector<int>& Nh0,
    std::function<TYPE(const Array<TYPE>&)> H,
    const std::vector<imfilter::BC>& BC0, const TYPE X )
{
    PROFILE_START( "imfilter (lambda)" );
    IMFILTER_ASSERT( A.ndim() == Nh0.size() );
    IMFILTER_ASSERT( A.ndim() == BC0.size() );
    std::vector<size_t> Nh2( A.size() );
    for (int d=0; d<A.ndim(); d++)
        Nh2[d] = 2*Nh0[d]+1;
    auto B = A;
    Array<TYPE> data(Nh2);
    IMFILTER_INSIST(A.ndim()<=3,"Not programmed for more than 3 dimensions yet");
    auto N = A.size();
    auto Nh = Nh0;
    auto BC = BC0;
    N.resize(3,1);
    Nh.resize(3,0);
    BC.resize(3,imfilter::BC::fixed);
    for ( int k1 = 0; k1 < N[2]; k1++ ) {
        for ( int j1 = 0; j1 < N[1]; j1++ ) {
            for ( int i1 = 0; i1 < N[0]; i1++ ) {
                for ( int kh = -Nh[2]; kh <= Nh[2]; kh++ ) {
                    int k2 = imfilter_index( k1+kh, N[2], BC[2] );
                    for ( int jh = -Nh[1]; jh <= Nh[1]; jh++ ) {
                        int j2 = imfilter_index( j1+jh, N[1], BC[1] );
                        for ( int ih = -Nh[0]; ih <= Nh[0]; ih++ ) {
                            int i2 = imfilter_index( i1+ih, N[0], BC[0] );
                            bool fixed = i2 == -1 || j2 == -1 || k2 == -1;
                            data(ih+Nh[0],jh+Nh[1],kh+Nh[2]) = fixed ? X : A(i2,j2,k2);
                        }
                    }
                }
                B(i1,j1,k1) = H( data );
            }
        }
    }
    PROFILE_STOP( "imfilter (lambda)" );
    return B;
}


/********************************************************
* imfilter with separable filter functions              *
********************************************************/
template<class TYPE>
Array<TYPE> imfilter::imfilter_separable( const Array<TYPE>& A,
    const std::vector<Array<TYPE>>& H,
    const std::vector<imfilter::BC>& boundary, const TYPE X )
{
    PROFILE_START( "imfilter_separable" );
    IMFILTER_ASSERT( A.ndim() == (int) H.size() );
    IMFILTER_ASSERT( A.ndim() == (int) boundary.size() );
    std::vector<size_t> Nh( H.size() );
    for (int d=0; d<A.ndim(); d++) {
        IMFILTER_ASSERT(H[d].ndim()==1);
        Nh[d] = (H[d].length()-1)/2;
        IMFILTER_INSIST(2*Nh[d]+1==H[d].length(),"Filter must be of size 2*N+1");
    }
    auto B = A;
    for ( int d = 0; d < A.ndim(); d++ ) {
        int N = A.size(d);
        int Ns = 1;
        int Ne = 1;
        for ( int d2 = 0; d2 < d; d2++ )
            Ns *= A.size(d2);
        for ( int d2 = d+1; d2 < A.ndim(); d2++ )
            Ne *= A.size(d2);
        filter_direction( Ns, N, Ne, Nh[d], H[d].data(), boundary[d], X, B.data() );
    }
    PROFILE_STOP( "imfilter_separable" );
    return B;
}
template<class TYPE>
Array<TYPE> imfilter::imfilter_separable( const Array<TYPE>& A, const std::vector<int>& Nh,
    std::vector<std::function<TYPE(const Array<TYPE>&)>> H,
    const std::vector<imfilter::BC>& boundary, const TYPE X )
{
    PROFILE_START( "imfilter_separable (lambda)" );
    IMFILTER_ASSERT( A.ndim() == (int) boundary.size() );
    auto B = A;
    for ( int d = 0; d < A.ndim(); d++ ) {
        int N = A.size(d);
        int Ns = 1;
        int Ne = 1;
        for ( int d2 = 0; d2 < d; d2++ )
            Ns *= A.size(d2);
        for ( int d2 = d+1; d2 < A.ndim(); d2++ )
            Ne *= A.size(d2);
        filter_direction( Ns, N, Ne, Nh[d], H[d], boundary[d], X, B.data() );
    }
    PROFILE_STOP( "imfilter_separable (lambda)" );
    return B;
}
template<class TYPE>
Array<TYPE> imfilter::imfilter_separable( const Array<TYPE>& A, const std::vector<int>& Nh,
    std::vector<std::function<TYPE(int, const TYPE*)>> H,
    const std::vector<imfilter::BC>& boundary, const TYPE X )
{
    PROFILE_START( "imfilter_separable (function)" );
    IMFILTER_ASSERT( A.ndim() == (int) boundary.size() );
    auto B = A;
    for ( int d = 0; d < A.ndim(); d++ ) {
        int N = A.size(d);
        int Ns = 1;
        int Ne = 1;
        for ( int d2 = 0; d2 < d; d2++ )
            Ns *= A.size(d2);
        for ( int d2 = d+1; d2 < A.ndim(); d2++ )
            Ne *= A.size(d2);
        filter_direction( Ns, N, Ne, Nh[d], H[d], boundary[d], X, B.data() );
    }
    PROFILE_STOP( "imfilter_separable (function)" );
    return B;
}



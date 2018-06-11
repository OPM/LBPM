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
#ifndef COMMUNICATION_HPP_INC
#define COMMUNICATION_HPP_INC

#include "common/Communication.h"
#include "common/MPI_Helpers.h"
#include "common/Utilities.h"
//#include "ProfilerApp.h"


/********************************************************
*  Structure to store the rank info                     *
********************************************************/
template<class TYPE>
fillHalo<TYPE>::fillHalo( MPI_Comm comm_, const RankInfoStruct& info_,
    std::array<int,3> n_, std::array<int,3> ng_, int tag0, int depth_,
    std::array<bool,3> fill, std::array<bool,3> periodic ):
    comm(comm_), info(info_), n(n_), ng(ng_), depth(depth_)
{
    if ( std::is_same<TYPE,double>() ) {
        N_type = 1;
        datatype = MPI_DOUBLE;
    } else if ( std::is_same<TYPE,float>() ) {
        N_type = 1;
        datatype = MPI_FLOAT;
    } else if ( sizeof(TYPE)%sizeof(double)==0 ) {
        N_type = sizeof(TYPE) / sizeof(double);
        datatype = MPI_DOUBLE;
    } else if ( sizeof(TYPE)%sizeof(float)==0 ) {
        N_type = sizeof(TYPE) / sizeof(float);
        datatype = MPI_FLOAT;
    } else {
        N_type = sizeof(TYPE);
        datatype = MPI_BYTE;
    }
    // Set the fill pattern
    memset(fill_pattern,0,sizeof(fill_pattern));
    if ( fill[0] ) {
        fill_pattern[0][1][1] = true;
        fill_pattern[2][1][1] = true;
        fill_pattern[1][0][1] = true;
        fill_pattern[1][2][1] = true;
        fill_pattern[1][1][0] = true;
        fill_pattern[1][1][2] = true;
    }
    if ( fill[1] ) {
        fill_pattern[0][0][1] = true;
        fill_pattern[0][2][1] = true;
        fill_pattern[2][0][1] = true;
        fill_pattern[2][2][1] = true;
        fill_pattern[0][1][0] = true;
        fill_pattern[0][1][2] = true;
        fill_pattern[2][1][0] = true;
        fill_pattern[2][1][2] = true;
        fill_pattern[1][0][0] = true;
        fill_pattern[1][0][2] = true;
        fill_pattern[1][2][0] = true;
        fill_pattern[1][2][2] = true;
    }
    if ( fill[2] ) {
        fill_pattern[0][0][0] = true;
        fill_pattern[0][0][2] = true;
        fill_pattern[0][2][0] = true;
        fill_pattern[0][2][2] = true;
        fill_pattern[2][0][0] = true;
        fill_pattern[2][0][2] = true;
        fill_pattern[2][2][0] = true;
        fill_pattern[2][2][2] = true;
    }
    // Remove communication for non-perioidic directions
    if ( !periodic[0] && info.ix==0 ) {
        for (int j=0; j<3; j++) {
            for (int k=0; k<3; k++)
                fill_pattern[0][j][k] = false;
        }
    }
    if ( !periodic[0] && info.ix==info.nx-1 ) {
        for (int j=0; j<3; j++) {
            for (int k=0; k<3; k++)
                fill_pattern[2][j][k] = false;
        }
    }
    if ( !periodic[1] && info.jy==0 ) {
        for (int i=0; i<3; i++) {
            for (int k=0; k<3; k++)
                fill_pattern[i][0][k] = false;
        }
    }
    if ( !periodic[1] && info.jy==info.ny-1 ) {
        for (int i=0; i<3; i++) {
            for (int k=0; k<3; k++)
                fill_pattern[i][2][k] = false;
        }
    }
    if ( !periodic[2] && info.kz==0 ) {
        for (int i=0; i<3; i++) {
            for (int j=0; j<3; j++)
                fill_pattern[i][j][0] = false;
        }
    }
    if ( !periodic[2] && info.kz==info.nz-1 ) {
        for (int i=0; i<3; i++) {
            for (int j=0; j<3; j++)
                fill_pattern[i][j][2] = false;
        }
    }
    // Determine the number of elements for each send/recv
    for (int i=0; i<3; i++) {
        int ni = (i-1)==0 ? n[0]:ng[0];
        for (int j=0; j<3; j++) {
            int nj = (j-1)==0 ? n[1]:ng[1];
            for (int k=0; k<3; k++) {
                int nk = (k-1)==0 ? n[2]:ng[2];
                if ( fill_pattern[i][j][k] )
                    N_send_recv[i][j][k] = ni*nj*nk;
                else
                    N_send_recv[i][j][k] = 0;
            }
        }
    }
    // Create send/recv buffers
    size_t N_mem=0;
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            for (int k=0; k<3; k++)
                N_mem += N_send_recv[i][j][k];
        }
    }
    mem = new TYPE[2*depth*N_mem];
    size_t index = 0;
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            for (int k=0; k<3; k++) {
                send[i][j][k] = &mem[index];
                index += depth*N_send_recv[i][j][k];
                recv[i][j][k] = &mem[index];
                index += depth*N_send_recv[i][j][k];
            }
        }
    }
    // Create the tags
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            for (int k=0; k<3; k++) {
                tag[i][j][k] = tag0 + i + j*3 + k*9;
            }
        }
    }

}
template<class TYPE>
fillHalo<TYPE>::~fillHalo( )
{
    delete [] mem;
}
template<class TYPE>
void fillHalo<TYPE>::fill( Array<TYPE>& data )
{
    //PROFILE_START("fillHalo::fill",1);
    int depth2 = data.size(3);
    ASSERT((int)data.size(0)==n[0]+2*ng[0]);
    ASSERT((int)data.size(1)==n[1]+2*ng[1]);
    ASSERT((int)data.size(2)==n[2]+2*ng[2]);
    ASSERT(depth2<=depth);
    ASSERT(data.ndim()==3||data.ndim()==4);
    // Start the recieves
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            for (int k=0; k<3; k++) {
                if ( !fill_pattern[i][j][k] )
                    continue;
                MPI_Irecv( recv[i][j][k], N_type*depth2*N_send_recv[i][j][k], datatype, 
                    info.rank[i][j][k], tag[2-i][2-j][2-k], comm, &recv_req[i][j][k] );
            }
        }
    }
    // Pack the src data and start the sends
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            for (int k=0; k<3; k++) {
                if ( !fill_pattern[i][j][k] )
                    continue;
                pack( data, i-1, j-1, k-1, send[i][j][k] );
                MPI_Isend( send[i][j][k], N_type*depth2*N_send_recv[i][j][k], datatype, 
                    info.rank[i][j][k], tag[i][j][k], comm, &send_req[i][j][k] );
            }
        }
    }
    // Recv the dst data and unpack (we recive in reverse order to match the sends)
    MPI_Status status;
    for (int i=2; i>=0; i--) {
        for (int j=2; j>=0; j--) {
            for (int k=2; k>=0; k--) {
                if ( !fill_pattern[i][j][k] )
                    continue;
                MPI_Wait(&recv_req[i][j][k],&status);
                unpack( data, i-1, j-1, k-1, recv[i][j][k] );
            }
        }
    }
    // Wait until all sends have completed
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            for (int k=0; k<3; k++) {
                if ( !fill_pattern[i][j][k] )
                    continue;
                MPI_Wait(&send_req[i][j][k],&status);
            }
        }
    }
    //PROFILE_STOP("fillHalo::fill",1);
}
template<class TYPE>
void fillHalo<TYPE>::pack( const Array<TYPE>& data, int i0, int j0, int k0, TYPE *buffer )
{
    int depth2 = data.size(3);
    int ni = i0==0 ? n[0]:ng[0];
    int nj = j0==0 ? n[1]:ng[1];
    int nk = k0==0 ? n[2]:ng[2];
    int is = i0==0 ? ng[0]:((i0==-1)?ng[0]:n[0]);
    int js = j0==0 ? ng[1]:((j0==-1)?ng[1]:n[1]);
    int ks = k0==0 ? ng[2]:((k0==-1)?ng[2]:n[2]);
    for (int d=0; d<depth2; d++) {
        for (int k=0; k<nk; k++) {
            for (int j=0; j<nj; j++) {
                for (int i=0; i<ni; i++) {
                    buffer[i+j*ni+k*ni*nj+d*ni*nj*nk] = data(i+is,j+js,k+ks,d);
                }
            }
        }
    }
}
template<class TYPE>
void fillHalo<TYPE>::unpack( Array<TYPE>& data, int i0, int j0, int k0, const TYPE *buffer )
{
    int depth2 = data.size(3);
    int ni = i0==0 ? n[0]:ng[0];
    int nj = j0==0 ? n[1]:ng[1];
    int nk = k0==0 ? n[2]:ng[2];
    int is = i0==0 ? ng[0]:((i0==-1)?0:n[0]+ng[0]);
    int js = j0==0 ? ng[1]:((j0==-1)?0:n[1]+ng[1]);
    int ks = k0==0 ? ng[2]:((k0==-1)?0:n[2]+ng[2]);
    for (int d=0; d<depth2; d++) {
        for (int k=0; k<nk; k++) {
            for (int j=0; j<nj; j++) {
                for (int i=0; i<ni; i++) {
                    data(i+is,j+js,k+ks,d) = buffer[i+j*ni+k*ni*nj+d*ni*nj*nk];
                }
            }
        }
    }
}


/********************************************************
*  Function to remove the ghost halo                    *
********************************************************/
template<class TYPE>
template<class TYPE1, class TYPE2>
void fillHalo<TYPE>::copy( const Array<TYPE1>& src, Array<TYPE2>& dst )
{
    //PROFILE_START("fillHalo::copy",1);
    ASSERT( (int)src.size(0)==n[0] || (int)src.size(0)==n[0]+2*ng[0] );
    ASSERT( (int)dst.size(0)==n[0] || (int)dst.size(0)==n[0]+2*ng[0] );
    bool src_halo = (int)src.size(0)==n[0]+2*ng[0];
    bool dst_halo = (int)dst.size(0)==n[0]+2*ng[0];
    if ( src_halo ) {
        ASSERT((int)src.size(0)==n[0]+2*ng[0]);
        ASSERT((int)src.size(1)==n[1]+2*ng[1]);
        ASSERT((int)src.size(2)==n[2]+2*ng[2]);
    } else {
        ASSERT((int)src.size(0)==n[0]);
        ASSERT((int)src.size(1)==n[1]);
        ASSERT((int)src.size(2)==n[2]);
    }
    if ( dst_halo ) {
        ASSERT((int)dst.size(0)==n[0]+2*ng[0]);
        ASSERT((int)dst.size(1)==n[1]+2*ng[1]);
        ASSERT((int)dst.size(2)==n[2]+2*ng[2]);
    } else {
        ASSERT((int)dst.size(0)==n[0]);
        ASSERT((int)dst.size(1)==n[1]);
        ASSERT((int)dst.size(2)==n[2]);
    }
    if ( src_halo == dst_halo ) {
        // Src and dst halos match
        for (size_t i=0; i<src.length(); i++)
            dst(i) = src(i);
    } else if ( src_halo && !dst_halo ) {
        // Src has halos
        for (int k=0; k<n[2]; k++) {
            for (int j=0; j<n[1]; j++) {
                for (int i=0; i<n[0]; i++) {
                    dst(i,j,k) = src(i+ng[0],j+ng[1],k+ng[2]);
                }
            }
        }
    } else if ( !src_halo && dst_halo ) {
        // Dst has halos
        for (int k=0; k<n[2]; k++) {
            for (int j=0; j<n[1]; j++) {
                for (int i=0; i<n[0]; i++) {
                    dst(i+ng[0],j+ng[1],k+ng[2]) = src(i,j,k);
                }
            }
        }
        fill(dst);
    }
    //PROFILE_STOP("fillHalo::copy",1);
}


#endif

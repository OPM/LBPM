#ifndef COMMUNICATION_HPP_INC
#define COMMUNICATION_HPP_INC

#include "common/Communication.h"
#include "common/MPI_Helpers.h"
#include "common/Utilities.h"
//#include "ProfilerApp.h"


/********************************************************
*  Redistribute data between two grids                  *
********************************************************/
template<class TYPE>
Array<TYPE> redistribute( const RankInfoStruct& src_rank, const Array<TYPE>& src_data,
    const RankInfoStruct& dst_rank, std::array<int,3> dst_size, MPI_Comm comm )
{
#ifdef USE_MPI
    // Get the src size
    std::array<int,3> src_size;
    int size0[3] = { (int) src_data.size(0), (int) src_data.size(1), (int) src_data.size(2) };
    MPI_Allreduce( size0, src_size.data(), 3, MPI_INT, MPI_MAX, comm );
    if ( !src_data.empty() )
        ASSERT( src_size[0] == size0[0] && src_size[1] == size0[1] && src_size[2] == size0[2] );
    // Check that dst_size matches on all ranks
    MPI_Allreduce( dst_size.data(), size0, 3, MPI_INT, MPI_MAX, comm );
    ASSERT( dst_size[0] == size0[0] && dst_size[1] == size0[1] && dst_size[2] == size0[2] );
    // Function to get overlap range
    auto calcOverlap = []( int i1[3], int i2[3], int j1[3], int j2[3] ) {
        std::vector<size_t> index;
        if ( i1[0] > j2[0] || i2[0] < j1[0] || i1[1] > j2[1] || i2[1] < j1[1] || i1[2] > j2[2] || i2[2] < j1[2] )
            return index;
        index.resize( 6 );
        index[0] = std::max( j1[0] - i1[0], 0 );
        index[1] = std::min( j2[0] - i1[0], i2[0] - i1[0] );
        index[2] = std::max( j1[1] - i1[1], 0 );
        index[3] = std::min( j2[1] - i1[1], i2[1] - i1[1] );
        index[4] = std::max( j1[2] - i1[2], 0 );
        index[5] = std::min( j2[2] - i1[2], i2[2] - i1[2] );
        return index;
    };
    // Pack and send my data to the appropriate ranks (including myself)
    std::vector<int> send_rank;
    std::vector<Array<TYPE>> send_data;
    if ( !src_data.empty() ) {
        int i1[3] = { src_size[0] * src_rank.ix, src_size[1] * src_rank.jy, src_size[2] * src_rank.kz };
        int i2[3] = { i1[0] + src_size[0] - 1, i1[1] + src_size[1] - 1, i1[2] + src_size[2] - 1 };
        for ( size_t i=0; i<dst_rank.nx; i++ ) {
            for ( size_t j=0; j<dst_rank.ny; j++ ) {
                for ( size_t k=0; k<dst_rank.nz; k++ ) {
                    int j1[3] = { i * dst_size[0], j * dst_size[1], k * dst_size[2] };
                    int j2[3] = { j1[0] + dst_size[0] - 1, j1[1] + dst_size[1] - 1, j1[2] + dst_size[2] - 1 };
                    auto index = calcOverlap( i1, i2, j1, j2 );
                    if ( index.empty() )
                        continue;
                    send_rank.push_back( dst_rank.getRankForBlock(i,j,k) );
                    send_data.push_back( src_data.subset( index ) );
                }
            }
        }
    }
    std::vector<MPI_Request> send_request( send_rank.size() );
    for (size_t i=0; i<send_rank.size(); i++)
        MPI_Isend( send_data[i].data(), sizeof(TYPE)*send_data[i].length(), MPI_BYTE, send_rank[i], 5462, comm, &send_request[i]);
    // Unpack data from the appropriate ranks (including myself)
    Array<TYPE> dst_data( dst_size[0], dst_size[1], dst_size[2] );
    int i1[3] = { dst_size[0] * dst_rank.ix, dst_size[1] * dst_rank.jy, dst_size[2] * dst_rank.kz };
    int i2[3] = { i1[0] + dst_size[0] - 1, i1[1] + dst_size[1] - 1, i1[2] + dst_size[2] - 1 };
    for ( size_t i=0; i<src_rank.nx; i++ ) {
        for ( size_t j=0; j<src_rank.ny; j++ ) {
            for ( size_t k=0; k<src_rank.nz; k++ ) {
                int j1[3] = { i * src_size[0], j * src_size[1], k * src_size[2] };
                int j2[3] = { j1[0] + src_size[0] - 1, j1[1] + src_size[1] - 1, j1[2] + src_size[2] - 1 };
                auto index = calcOverlap( i1, i2, j1, j2 );
                if ( index.empty() )
                    continue;
                int rank  = src_rank.getRankForBlock(i,j,k);
                Array<TYPE> data( index[1] - index[0] + 1, index[3] - index[2] + 1, index[5] - index[4] + 1 );
                MPI_Recv( data.data(), sizeof(TYPE)*data.length(), MPI_BYTE, rank, 5462, comm, MPI_STATUS_IGNORE );
                dst_data.copySubset( index, data );
            }
        }
    }
    // Free data
    MPI_Waitall( send_request.size(), send_request.data(), MPI_STATUSES_IGNORE );
    return dst_data;
#else
    return src_data.subset( { 0, dst_size[0]-1, 0, dst_size[1]-1, 0, dst_size[2]-1 );
#endif
}



/********************************************************
*  Structure to fill halo cells                         *
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

/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University
  Copyright Equnior ASA

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
#ifndef Domain_INC
#define Domain_INC

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <time.h>
#include <exception>
#include <stdexcept>

#include "common/Array.h"
#include "common/Utilities.h"
#include "common/MPI.h"
#include "common/Communication.h"
#include "common/Database.h"
/**
 * @file Domain.h
 * \brief Parallel Domain data structures and helper functions
 */

/**
 * \class Box
 *
 * @details
 * information about a box
 */
class Box {
public:
    int ifirst[3];
    int ilast[3];
};

class Patch;

/**
 * \class Domain
 *
 * @details
 * the Domain class includes basic information to distribution 3D image data to multiple processes using MPI.
 * A regular domain decomposision is performed, with each MPI process getting a [Nx,Ny,Nz] sub-domain.
 * 8-bit image labels are retained internally.
 * The domain class resides on the CPU and provides utilities to support CPU-based analysis.
 * GPU-based data structures should be constructed separately but may utilize information that the Domain class provides.
*/

class Domain {
public:
    /**
    * \brief Constructor
    * @param db           input database
    * @param Communicator MPI communicator  
    */
    Domain(std::shared_ptr<Database> db, const Utilities::MPI &Communicator);

    /**
     * \brief Obsolete constructor
     */
    Domain(int nx, int ny, int nz, int rnk, int npx, int npy, int npz,
           double lx, double ly, double lz, int BC);

    /**
    * \brief Empty constructor
    */
    Domain() = delete;

    /**
    * \brief Copy constructor
    */
    Domain(const Domain &) = delete;

    /**
     * \brief Assignment operator
     */
    Domain &operator=(const Domain &) = delete;

    /**
     * \brief Destructor
     */
    ~Domain();

    /**
     * \brief Get the database
    */
    inline std::shared_ptr<const Database> getDatabase() const { return d_db; }

    /** 
    * \brief Get the domain box
    */
    inline const Box &getBox() const { return d_box; }

    /** 
    * \brief Get local patch
    */
    inline const Patch &getLocalPatch() const { return *d_localPatch; }

    /** 
     * \brief Get all patches
    */
    inline const std::vector<Patch> &getAllPatch() const { return d_patches; }

private:
    /** 
     * \brief initialize from database
    */
    void initialize(std::shared_ptr<Database> db);

    std::shared_ptr<Database> d_db;
    Box d_box;
    Patch *d_localPatch;
    std::vector<Patch> d_patches;

public: // Public variables (need to create accessors instead)
    std::shared_ptr<Database> database;
    double Lx, Ly, Lz, Volume, voxel_length;
    int Nx, Ny, Nz, N;
    int inlet_layers_x, inlet_layers_y, inlet_layers_z;
    int outlet_layers_x, outlet_layers_y, outlet_layers_z;
    int offset_x, offset_y, offset_z;
    int inlet_layers_phase; //as usual: 1->n, 2->w
    int outlet_layers_phase;
    double porosity;
    RankInfoStruct rank_info;

    Utilities::MPI Comm; // MPI Communicator for this domain

    int BoundaryCondition;

    //**********************************
    // MPI ranks for all 18 neighbors
    //**********************************
    /** 
     * \brief Compute the porosity based on the current domain id file
    */
    inline double Porosity() const { return porosity; }
    inline int iproc() const { return rank_info.ix; }
    inline int jproc() const { return rank_info.jy; }
    inline int kproc() const { return rank_info.kz; }
    inline int nprocx() const { return rank_info.nx; }
    inline int nprocy() const { return rank_info.ny; }
    inline int nprocz() const { return rank_info.nz; }
    inline int rank() const { return rank_info.rank[1][1][1]; }
    inline int rank_X() const { return rank_info.rank[2][1][1]; }
    inline int rank_x() const { return rank_info.rank[0][1][1]; }
    inline int rank_Y() const { return rank_info.rank[1][2][1]; }
    inline int rank_y() const { return rank_info.rank[1][0][1]; }
    inline int rank_Z() const { return rank_info.rank[1][1][2]; }
    inline int rank_z() const { return rank_info.rank[1][1][0]; }
    inline int rank_XY() const { return rank_info.rank[2][2][1]; }
    inline int rank_xy() const { return rank_info.rank[0][0][1]; }
    inline int rank_Xy() const { return rank_info.rank[2][0][1]; }
    inline int rank_xY() const { return rank_info.rank[0][2][1]; }
    inline int rank_XZ() const { return rank_info.rank[2][1][2]; }
    inline int rank_xz() const { return rank_info.rank[0][1][0]; }
    inline int rank_Xz() const { return rank_info.rank[2][1][0]; }
    inline int rank_xZ() const { return rank_info.rank[0][1][2]; }
    inline int rank_YZ() const { return rank_info.rank[1][2][2]; }
    inline int rank_yz() const { return rank_info.rank[1][0][0]; }
    inline int rank_Yz() const { return rank_info.rank[1][2][0]; }
    inline int rank_yZ() const { return rank_info.rank[1][0][2]; }

    //**********************************
    //......................................................................................
    // Get the actual D3Q19 communication counts (based on location of solid phase)
    // Discrete velocity set symmetry implies the sendcount = recvcount
    //......................................................................................
    inline int recvCount(const char *dir) const {
        return getRecvList(dir).size();
    }
    inline int sendCount(const char *dir) const {
        return getSendList(dir).size();
    }
    inline const int *recvList(const char *dir) const {
        return getRecvList(dir).data();
    }
    inline const int *sendList(const char *dir) const {
        return getSendList(dir).data();
    }

    //......................................................................................
    // Solid indicator function
    std::vector<signed char> id;

    /** 
     * \brief Read domain IDs from file
    */
    void ReadIDs();

    /** 
     * \brief Read domain IDs from SWC file
    */
    void read_swc(const std::string &Filename);

    /** 
     * \brief Compute the porosity
    */
    void ComputePorosity();

    /** 
     * \brief Read image and perform domain decomposition
     * @param filename  - name of file to read IDs
    */
    void Decomp(const std::string &filename);

    /** 
     * \brief Perform a halo exchange using MPI
     * @param Mesh - array data that holds scalar values
    */
    void CommunicateMeshHalo(DoubleArray &Mesh);

    /** 
     * \brief Initialize communication data structures within Domain object. 
     * This routine needs to be called before the communication functionality will work
    */
    void CommInit();

    /** 
     * \brief Count number of pore nodes (labels > 1)
    */
    int PoreCount();

    /** 
     * \brief Read array data from a file and distribute to local arrays for each MPI process
     * @param Filename - name of the file to read the data 
     * @param Datatype - data type to use 
     * @param UserData - Array to store the values that are read
    */
    void ReadFromFile(const std::string &Filename, const std::string &Datatype,
                      double *UserData);

    /**
     * \brief Aggregate labels from all MPI processes and write to a file
     * @param filename - name of the file to write
     */
    void AggregateLabels(const std::string &filename);
    /**
     * \brief Aggregate user provided array  from all MPI processes and write to a single file
     * @param filename - name of the file to write
     * @param UserData - array data to aggregate and write
     */
    void AggregateLabels(const std::string &filename, DoubleArray &UserData);

private:
    /**
     * \brief Pack halo data for 8-bit integer
     * @param list - list of values in the halo
     * @param count - count of values in the halo
     * @param sendbuf - memory buffer to use to pack values for MPI
     * @param ID - 8-bit values on mesh [Nx, Ny, Nz]
     */
    void PackID(int *list, int count, signed char *sendbuf, signed char *ID);

    /**
     * \brief Unpack halo data for 8-bit integer
     * @param list - list of values in the halo
     * @param count - count of values in the halo
     * @param recvbuf - memory buffer containing values recieved by MPI
     * @param ID - 8-bit values on mesh [Nx, Ny, Nz]
     */
    void UnpackID(int *list, int count, signed char *recvbuf, signed char *ID);

    //......................................................................................
    MPI_Request req1[18], req2[18];
    //......................................................................................
    std::vector<int> sendList_x, sendList_y, sendList_z, sendList_X, sendList_Y,
        sendList_Z;
    std::vector<int> sendList_xy, sendList_yz, sendList_xz, sendList_Xy,
        sendList_Yz, sendList_xZ;
    std::vector<int> sendList_xY, sendList_yZ, sendList_Xz, sendList_XY,
        sendList_YZ, sendList_XZ;
    //......................................................................................
    std::vector<int> recvList_x, recvList_y, recvList_z, recvList_X, recvList_Y,
        recvList_Z;
    std::vector<int> recvList_xy, recvList_yz, recvList_xz, recvList_Xy,
        recvList_Yz, recvList_xZ;
    std::vector<int> recvList_xY, recvList_yZ, recvList_Xz, recvList_XY,
        recvList_YZ, recvList_XZ;
    //......................................................................................
    const std::vector<int> &getRecvList(const char *dir) const;
    const std::vector<int> &getSendList(const char *dir) const;
};

template <class TYPE> class PatchData;

enum class DataLocation { CPU, DEVICE };

/**
 * \class Patch
 *
 * @details
 * store patch data
 */
class Patch {
public:
    //! Empty constructor
    Patch() = delete;

    //! Copy constructor
    Patch(const Patch &) = delete;

    //! Assignment operator
    Patch &operator=(const Patch &) = delete;

    //! Return the box for the patch
    inline const Box &getBox() const { return d_box; }

    //! Create patch data
    template <class TYPE>
    std::shared_ptr<PatchData<TYPE>>
    createPatchData(DataLocation location) const;

private:
    Box d_box;
    int d_owner;
    Domain *d_domain;
};

// Class to hold data on a patch
template <class TYPE> class PatchData {
public:
    //! Get the raw data pointer
    TYPE *data() { return d_data; }

    //! Get the raw data pointer
    const TYPE *data() const { return d_data; }

    //! Get the patch
    const Patch &getPatch() const { return *d_patch; }

    //! Start communication
    void beginCommunication();

    //! End communication
    void endCommunication();

    //! Access ghost values
    TYPE operator()(int, int, int) const;

    //! Copy data from another PatchData
    void copy(const PatchData &rhs);

private:
    DataLocation d_location;
    const Patch *d_patch;
    TYPE *d_data;
    TYPE *d_gcw;
};

void WriteCheckpoint(const char *FILENAME, const double *cDen,
                     const double *cfq, size_t Np);

void ReadCheckpoint(char *FILENAME, double *cDen, double *cfq, size_t Np);

void ReadBinaryFile(char *FILENAME, double *Data, size_t N);

#endif

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
// This file impliments a wrapper class for MPI functions

#include "common/MPI.h"
#include "common/Utilities.h"
#include "common/Utilities.hpp"

#include "ProfilerApp.h"
#include "StackTrace/ErrorHandlers.h"
#include "StackTrace/StackTrace.h"

// Include all other headers
#include <algorithm>
#include <chrono>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <random>
#include <stdexcept>
#include <thread>
#include <typeinfo>

// Include OS specific headers
#undef USE_WINDOWS
#undef USE_LINUX
#undef USE_MAC
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
// We are using windows
#define USE_WINDOWS
#include <process.h>
#include <windows.h>
#define sched_yield() Sleep(0)
#elif defined(__APPLE__)
// Using MAC
#define USE_MAC
#include <sched.h>
#elif defined(__linux) || defined(__linux__) || defined(__unix) ||             \
    defined(__posix)
// We are using linux
#define USE_LINUX
#include <sched.h>
#include <unistd.h>
#else
#error Unknown OS
#endif

// Convience defines
#define MPI_ERROR ERROR
#define MPI_ASSERT ASSERT
#define MPI_INSIST INSIST
#define MPI_WARNING WARNING
#define MPI_CLASS_COMM_NULL MPI_COMM_NULL
#define MPI_CLASS_COMM_SELF MPI_COMM_SELF
#define MPI_CLASS_COMM_WORLD MPI_COMM_WORLD

// Global variable to track create new unique comms (dup and split)
#ifndef USE_MPI
MPI_Comm uniqueGlobalComm = 11;
#endif

#if defined(USE_SAMRAI) && defined(USE_PETSC) && !defined(USE_MPI)
int MPI_REQUEST_NULL = 3;
int MPI_ERR_IN_STATUS = 4;
#endif

namespace Utilities {

// Some special structs to work with MPI
#ifdef USE_MPI
struct IntIntStruct {
    int j;
    int i;
};
struct LongIntStruct {
    long int j;
    int i;
};
struct FloatIntStruct {
    float f;
    int i;
};
struct DoubleIntStruct {
    double d;
    int i;
};
#endif

// Initialized the static member variables
volatile unsigned int MPI_CLASS::N_MPI_Comm_created = 0;
volatile unsigned int MPI_CLASS::N_MPI_Comm_destroyed = 0;
short MPI_CLASS::profile_level = 127;

// Define a type for use with size_t
#ifdef USE_MPI
static MPI_Datatype MPI_SIZE_T = 0x0;
static MPI_Datatype getSizeTDataType() {
    int size_int, size_long, size_longlong, size_longlong2;
    MPI_Type_size(MPI_UNSIGNED, &size_int);
    MPI_Type_size(MPI_UNSIGNED_LONG, &size_long);
    MPI_Type_size(MPI_UNSIGNED_LONG_LONG, &size_longlong);
    MPI_Type_size(MPI_LONG_LONG_INT, &size_longlong2);
    if (sizeof(size_t) == size_int) {
        return MPI_UNSIGNED;
    } else if (sizeof(size_t) == size_long) {
        return MPI_UNSIGNED_LONG;
    } else if (sizeof(size_t) == size_longlong) {
        return MPI_UNSIGNED_LONG_LONG;
    } else if (sizeof(size_t) == size_longlong2) {
        MPI_WARNING("Using signed long long datatype for size_t in MPI");
        return MPI_LONG_LONG_INT; // Note: this is not unsigned
    } else {
        MPI_ERROR("No suitable datatype found");
    }
    return 0;
}
#endif

// Static data for asyncronous communication without MPI
// Note: these routines may not be thread-safe yet
#ifndef USE_MPI
static const int mpi_max_tag = 0x003FFFFF;
struct Isendrecv_struct {
    const char *data; // Pointer to data
    int status;       // Status: 1-sending, 2-recieving
};
std::map<MPI_Request, Isendrecv_struct> global_isendrecv_list;
static MPI_Request getRequest(MPI_Comm comm, int tag) {
    MPI_ASSERT(tag >= 0 && tag <= mpi_max_tag);
    // Use hashing function: 2^64*0.5*(sqrt(5)-1)
    uint64_t a = static_cast<uint8_t>(comm) * 0x9E3779B97F4A7C15;
    uint64_t b = static_cast<uint8_t>(tag) * 0x9E3779B97F4A7C15;
    uint64_t hash = a ^ b;
    MPI_Request request;
    memcpy(&request, &hash, sizeof(MPI_Request));
    return request;
}
#endif

// Check the mpi error code
#ifdef USE_MPI
inline void check_MPI(int error) {
    if (error != MPI_SUCCESS)
        MPI_ERROR("Error calling MPI routine");
}
#endif

/******************************************************************
 * Some helper functions to convert between signed/unsigned types  *
 ******************************************************************/
DISABLE_WARNINGS
static inline constexpr unsigned int offset_int() {
    return ~static_cast<unsigned int>(std::numeric_limits<int>::min()) + 1;
}
static inline constexpr unsigned long int offset_long() {
    return ~static_cast<long int>(std::numeric_limits<long int>::min()) + 1;
}
static inline constexpr unsigned long long int offset_long_long() {
    return ~static_cast<long long int>(
               std::numeric_limits<long long int>::min()) +
           1;
}
ENABLE_WARNINGS
static inline unsigned int signed_to_unsigned(int x) {
    const auto offset = offset_int();
    return (x >= 0) ? static_cast<unsigned int>(x) + offset
                    : offset - static_cast<unsigned int>(-x);
}
static inline unsigned long int signed_to_unsigned(long int x) {
    const auto offset = offset_long();
    return (x >= 0) ? static_cast<unsigned long int>(x) + offset
                    : offset - static_cast<unsigned long int>(-x);
}
static inline unsigned long long int signed_to_unsigned(long long int x) {
    const auto offset = offset_long_long();
    return (x >= 0) ? static_cast<unsigned long long int>(x) + offset
                    : offset - static_cast<unsigned long long int>(-x);
}
static inline int unsigned_to_signed(unsigned int x) {
    const auto offset = offset_int();
    return (x >= offset) ? static_cast<int>(x - offset)
                         : -static_cast<int>(offset - x);
}
static inline long int unsigned_to_signed(unsigned long int x) {
    const auto offset = offset_long();
    return (x >= offset) ? static_cast<long int>(x - offset)
                         : -static_cast<long int>(offset - x);
}
static inline long long int unsigned_to_signed(unsigned long long int x) {
    const auto offset = offset_long_long();
    return (x >= offset) ? static_cast<long long int>(x - offset)
                         : -static_cast<long long int>(offset - x);
}

/************************************************************************
 *  Get the MPI version                                                  *
 ************************************************************************/
std::array<int, 2> MPI_CLASS::version() {
#ifdef USE_MPI
    int MPI_version;
    int MPI_subversion;
    MPI_Get_version(&MPI_version, &MPI_subversion);
    return {MPI_version, MPI_subversion};
#else
    return {0, 0};
#endif
}
std::string MPI_CLASS::info() {
#ifdef USE_MPI
#if MPI_VERSION >= 3
    int MPI_version_length = 0;
    char MPI_version_string[MPI_MAX_LIBRARY_VERSION_STRING];
    MPI_Get_library_version(MPI_version_string, &MPI_version_length);
    if (MPI_version_length > 0) {
        std::string MPI_info(MPI_version_string, MPI_version_length);
        size_t pos = MPI_info.find('\n');
        while (pos != std::string::npos) {
            MPI_info.insert(pos + 1, "   ");
            pos = MPI_info.find('\n', pos + 1);
        }
        return MPI_info;
    }
#endif
    auto tmp = version();
    return std::to_string(tmp[0]) + "." + std::to_string(tmp[0]);
#else
    return std::string();
#endif
}

/************************************************************************
 *  Functions to get/set the process affinities                          *
 ************************************************************************/
int MPI_CLASS::getNumberOfProcessors() {
    return std::thread::hardware_concurrency();
}
std::vector<int> MPI_CLASS::getProcessAffinity() {
    std::vector<int> procs;
#ifdef USE_LINUX
    cpu_set_t mask;
    int error = sched_getaffinity(getpid(), sizeof(cpu_set_t), &mask);
    if (error != 0)
        MPI_ERROR("Error getting process affinity");
    for (int i = 0; i < (int)sizeof(cpu_set_t) * CHAR_BIT; i++) {
        if (CPU_ISSET(i, &mask))
            procs.push_back(i);
    }
#elif defined(USE_MAC)
    // MAC does not support getting or setting the affinity
    printf("Warning: MAC does not support getting the process affinity\n");
    procs.clear();
#elif defined(USE_WINDOWS)
    HANDLE hProc = GetCurrentProcess();
    size_t procMask;
    size_t sysMask;
    PDWORD_PTR procMaskPtr = reinterpret_cast<PDWORD_PTR>(&procMask);
    PDWORD_PTR sysMaskPtr = reinterpret_cast<PDWORD_PTR>(&sysMask);
    GetProcessAffinityMask(hProc, procMaskPtr, sysMaskPtr);
    for (int i = 0; i < (int)sizeof(size_t) * CHAR_BIT; i++) {
        if ((procMask & 0x1) != 0)
            procs.push_back(i);
        procMask >>= 1;
    }
#else
#error Unknown OS
#endif
    return procs;
}
void MPI_CLASS::setProcessAffinity(const std::vector<int> &procs) {
#ifdef USE_LINUX
    cpu_set_t mask;
    CPU_ZERO(&mask);
    for (auto cpu : procs)
        CPU_SET(cpu, &mask);
    int error = sched_setaffinity(getpid(), sizeof(cpu_set_t), &mask);
    if (error != 0)
        MPI_ERROR("Error setting process affinity");
#elif defined(USE_MAC)
    // MAC does not support getting or setting the affinity
    NULL_USE(procs);
#elif defined(USE_WINDOWS)
    DWORD mask = 0;
    for (size_t i = 0; i < procs.size(); i++)
        mask |= ((DWORD)1) << procs[i];
    HANDLE hProc = GetCurrentProcess();
    SetProcessAffinityMask(hProc, mask);
#else
#error Unknown OS
#endif
}

/************************************************************************
 *  Function to check if MPI is active                                   *
 ************************************************************************/
bool MPI_CLASS::MPI_active() {
#ifdef USE_MPI
    int initialized = 0, finalized = 0;
    MPI_Initialized(&initialized);
    MPI_Finalized(&finalized);
    return initialized != 0 && finalized == 0;
#else
    return true;
#endif
}
MPI_CLASS::ThreadSupport MPI_CLASS::queryThreadSupport() {
#ifdef USE_MPI
    int provided = 0;
    MPI_Query_thread(&provided);
    if (provided == MPI_THREAD_SINGLE)
        return ThreadSupport::SINGLE;
    if (provided == MPI_THREAD_FUNNELED)
        return ThreadSupport::FUNNELED;
    if (provided == MPI_THREAD_SERIALIZED)
        return ThreadSupport::SERIALIZED;
    if (provided == MPI_THREAD_MULTIPLE)
        return ThreadSupport::MULTIPLE;
    return ThreadSupport::SINGLE;
#else
    return ThreadSupport::MULTIPLE;
#endif
}

/************************************************************************
 *  Function to perform a load balance of the given processes            *
 ************************************************************************/
void MPI_CLASS::balanceProcesses(const MPI_CLASS &globalComm, const int method,
                                 const std::vector<int> &procs,
                                 const int N_min_in, const int N_max_in) {
    // Build the list of processors to use
    std::vector<int> cpus = procs;
    if (cpus.empty()) {
        for (int i = 0; i < getNumberOfProcessors(); i++)
            cpus.push_back(i);
    }
    // Handle the "easy cases"
    if (method == 1) {
        // Trivial case where we do not need any communication
        setProcessAffinity(cpus);
        return;
    }
    // Get the sub-communicator for the current node
    MPI_CLASS nodeComm = globalComm.splitByNode();
    int N_min = std::min<int>(std::max<int>(N_min_in, 1), cpus.size());
    int N_max = N_max_in;
    if (N_max == -1)
        N_max = cpus.size();
    N_max = std::min<int>(N_max, cpus.size());
    MPI_ASSERT(N_max >= N_min);
    // Perform the load balance within the node
    if (method == 2) {
        int N_proc = cpus.size() / nodeComm.getSize();
        N_proc = std::max<int>(N_proc, N_min);
        N_proc = std::min<int>(N_proc, N_max);
        std::vector<int> cpus2(N_proc, -1);
        for (int i = 0; i < N_proc; i++)
            cpus2[i] = cpus[(nodeComm.getRank() * N_proc + i) % cpus.size()];
        setProcessAffinity(cpus2);
    } else {
        MPI_ERROR("Unknown method for load balance");
    }
}

/************************************************************************
 *  Empty constructor                                                    *
 ************************************************************************/
MPI_CLASS::MPI_CLASS() {
// Initialize the data members to a defaul communicator of self
#ifdef USE_MPI
    communicator = MPI_COMM_NULL;
    d_maxTag = 0x7FFFFFFF;
#else
    communicator = MPI_CLASS_COMM_NULL;
    d_maxTag = mpi_max_tag;
#endif
    d_count = nullptr;
    d_manage = false;
    comm_rank = 0;
    comm_size = 1;
    d_isNull = true;
    d_currentTag = nullptr;
    d_call_abort = true;
    tmp_alignment = -1;
}

/************************************************************************
 *  Empty deconstructor                                                  *
 ************************************************************************/
MPI_CLASS::~MPI_CLASS() { reset(); }
void MPI_CLASS::reset() {
    // Decrement the count if used
    int count = -1;
    if (d_count != nullptr)
        count = --(*d_count);
    if (count == 0) {
        // We are holding that last reference to the MPI_Comm object, we need to free it
        if (d_manage) {
#ifdef USE_MPI
            MPI_Comm_set_errhandler(communicator, MPI_ERRORS_ARE_FATAL);
            int err = MPI_Comm_free(&communicator);
            if (err != MPI_SUCCESS)
                MPI_ERROR("Problem free'ing MPI_Comm object");
            communicator = MPI_CLASS_COMM_NULL;
            ++N_MPI_Comm_destroyed;
#endif
        }
        delete d_count;
    }
    if (d_currentTag == nullptr) {
        // No tag index
    } else if (d_currentTag[1] > 1) {
        --(d_currentTag[1]);
    } else {
        delete[] d_currentTag;
    }
    d_manage = false;
    d_count = nullptr;
    comm_rank = 0;
    comm_size = 1;
    d_maxTag = 0;
    d_isNull = true;
    d_currentTag = nullptr;
    d_call_abort = true;
}

/************************************************************************
 *  Copy constructors                                                    *
 ************************************************************************/
MPI_CLASS::MPI_CLASS(const MPI_CLASS &comm)
    : communicator(comm.communicator), d_isNull(comm.d_isNull),
      d_manage(comm.d_manage), comm_rank(comm.comm_rank),
      comm_size(comm.comm_size), d_maxTag(comm.d_maxTag),
      d_currentTag(comm.d_currentTag) {
    // Initialize the data members to the existing comm object
    if (d_currentTag != nullptr)
        ++d_currentTag[1];
    d_call_abort = comm.d_call_abort;
    // Set and increment the count
    d_count = comm.d_count;
    if (d_count != nullptr)
        ++(*d_count);
    tmp_alignment = -1;
}
MPI_CLASS::MPI_CLASS(MPI_CLASS &&rhs) : MPI_CLASS() {
    std::swap(communicator, rhs.communicator);
    std::swap(d_isNull, rhs.d_isNull);
    std::swap(d_manage, rhs.d_manage);
    std::swap(d_call_abort, rhs.d_call_abort);
    std::swap(profile_level, rhs.profile_level);
    std::swap(comm_rank, rhs.comm_rank);
    std::swap(comm_size, rhs.comm_size);
    std::swap(d_maxTag, rhs.d_maxTag);
    std::swap(d_currentTag, rhs.d_currentTag);
    std::swap(d_count, rhs.d_count);
    std::swap(tmp_alignment, rhs.tmp_alignment);
}

/************************************************************************
 *  Assignment operators                                                 *
 ************************************************************************/
MPI_CLASS &MPI_CLASS::operator=(const MPI_CLASS &comm) {
    if (this == &comm) // protect against invalid self-assignment
        return *this;
    // Destroy the previous object
    this->reset();
    // Initialize the data members to the existing object
    this->communicator = comm.communicator;
    this->comm_rank = comm.comm_rank;
    this->comm_size = comm.comm_size;
    this->d_isNull = comm.d_isNull;
    this->d_manage = comm.d_manage;
    this->d_maxTag = comm.d_maxTag;
    this->d_call_abort = comm.d_call_abort;
    this->d_currentTag = comm.d_currentTag;
    if (this->d_currentTag != nullptr)
        ++(this->d_currentTag[1]);
    // Set and increment the count
    this->d_count = comm.d_count;
    if (this->d_count != nullptr)
        ++(*d_count);
    this->tmp_alignment = -1;
    return *this;
}
MPI_CLASS &MPI_CLASS::operator=(MPI_CLASS &&rhs) {
    if (this == &rhs) // protect against invalid self-assignment
        return *this;
    std::swap(communicator, rhs.communicator);
    std::swap(d_isNull, rhs.d_isNull);
    std::swap(d_manage, rhs.d_manage);
    std::swap(d_call_abort, rhs.d_call_abort);
    std::swap(profile_level, rhs.profile_level);
    std::swap(comm_rank, rhs.comm_rank);
    std::swap(comm_size, rhs.comm_size);
    std::swap(d_maxTag, rhs.d_maxTag);
    std::swap(d_currentTag, rhs.d_currentTag);
    std::swap(d_count, rhs.d_count);
    std::swap(tmp_alignment, rhs.tmp_alignment);
    return *this;
}

/************************************************************************
 *  Constructor from existing MPI communicator                           *
 ************************************************************************/
int d_global_currentTag_world1[2] = {1, 1};
int d_global_currentTag_world2[2] = {1, 1};
int d_global_currentTag_self[2] = {1, 1};
#ifdef USE_MPI
std::atomic_int d_global_count_world1 = {1};
std::atomic_int d_global_count_world2 = {1};
std::atomic_int d_global_count_self = {1};
#endif
MPI_CLASS::MPI_CLASS(MPI_Comm comm, bool manage) {
    d_count = nullptr;
    d_manage = false;
    tmp_alignment = -1;
    // Check if we are using our version of comm_world
    if (comm == MPI_CLASS_COMM_WORLD) {
        communicator = MPI_COMM_WORLD;
    } else if (comm == MPI_CLASS_COMM_SELF) {
        communicator = MPI_COMM_SELF;
    } else if (comm == MPI_CLASS_COMM_NULL) {
        communicator = MPI_COMM_NULL;
    } else {
        communicator = comm;
    }
#ifdef USE_MPI
    // We are using MPI, use the MPI communicator to initialize the data
    if (communicator != MPI_COMM_NULL) {
        // Set the MPI_SIZE_T datatype if it has not been set
        if (MPI_SIZE_T == 0x0)
            MPI_SIZE_T = getSizeTDataType();
        // Attach the error handler
        StackTrace::setMPIErrorHandler(communicator);
        // Get the communicator properties
        MPI_Comm_rank(communicator, &comm_rank);
        MPI_Comm_size(communicator, &comm_size);
        int flag, *val;
        int ierr = MPI_Comm_get_attr(communicator, MPI_TAG_UB, &val, &flag);
        MPI_ASSERT(ierr == MPI_SUCCESS);
        if (flag == 0) {
            d_maxTag =
                0x7FFFFFFF; // The tag is not a valid attribute (set to 2^31-1)
        } else {
            d_maxTag = *val;
            if (d_maxTag < 0) {
                d_maxTag = 0x7FFFFFFF;
            } // The maximum tag is > a signed int (set to 2^31-1)
            MPI_INSIST(d_maxTag >= 0x7FFF,
                       "maximum tag size is < MPI standard");
        }
    } else {
        comm_rank = 1;
        comm_size = 0;
        d_maxTag = 0x7FFFFFFF;
    }
    d_isNull = communicator == MPI_COMM_NULL;
    if (manage && communicator != MPI_COMM_NULL &&
        communicator != MPI_COMM_SELF && communicator != MPI_COMM_WORLD)
        d_manage = true;
    // Create the count (Note: we do not need to worry about thread safety)
    if (communicator == MPI_CLASS_COMM_WORLD) {
        d_count = &d_global_count_world1;
        ++(*d_count);
    } else if (communicator == MPI_COMM_WORLD) {
        d_count = &d_global_count_world2;
        ++(*d_count);
    } else if (communicator == MPI_COMM_SELF) {
        d_count = &d_global_count_self;
        ++(*d_count);
    } else if (communicator == MPI_COMM_NULL) {
        d_count = nullptr;
    } else {
        d_count = new std::atomic_int;
        *d_count = 1;
    }
    if (d_manage)
        ++N_MPI_Comm_created;

#else
    // We are not using MPI, intialize based on the communicator
    NULL_USE(manage);
    comm_rank = 0;
    comm_size = 1;
    d_maxTag = mpi_max_tag;
    d_isNull = communicator == MPI_COMM_NULL;
    if (d_isNull)
        comm_size = 0;
#endif
    if (communicator == MPI_CLASS_COMM_WORLD) {
        d_currentTag = d_global_currentTag_world1;
        ++(this->d_currentTag[1]);
    } else if (communicator == MPI_COMM_WORLD) {
        d_currentTag = d_global_currentTag_world2;
        ++(this->d_currentTag[1]);
    } else if (communicator == MPI_COMM_SELF) {
        d_currentTag = d_global_currentTag_self;
        ++(this->d_currentTag[1]);
    } else if (communicator == MPI_COMM_NULL) {
        d_currentTag = nullptr;
    } else {
        d_currentTag = new int[2];
        d_currentTag[0] = (d_maxTag <= 0x10000) ? 1 : 0x1FFF;
        d_currentTag[1] = 1;
    }
    d_call_abort = true;
}

/************************************************************************
 *  Return the ranks of the communicator in the global comm              *
 ************************************************************************/
std::vector<int> MPI_CLASS::globalRanks() const {
    if (d_isNull)
        return std::vector<int>();
#ifdef USE_MPI
    // Get my global rank and size if it has not been set
    static int globalRank = -1;
    static int globalSize = -1;
    if (globalRank == -1 && MPI_active()) {
        MPI_Comm_rank(MPI_CLASS_COMM_WORLD, &globalRank);
        MPI_Comm_size(MPI_CLASS_COMM_WORLD, &globalSize);
    }
    // Check if we are dealing with a serial or global communicator
    if (comm_size == 1)
        return std::vector<int>(1, globalRank);
    if (comm_size == globalSize) {
        std::vector<int> ranks(globalSize);
        for (int i = 0; i < globalSize; i++)
            ranks[i] = i;
        return ranks;
    }
    // Get the global rank from each rank in the communicator
    auto ranks = allGather(globalRank);
    std::sort(ranks.begin(), ranks.end());
    return ranks;
#else
    return std::vector<int>(1, 1);
#endif
}

/************************************************************************
 *  Generate a random number                                             *
 ************************************************************************/
size_t MPI_CLASS::rand() const {
    size_t val = 0;
    if (getRank() == 0) {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_int_distribution<size_t> dist;
        val = dist(gen);
    }
    val = bcast(val, 0);
    return val;
}

/************************************************************************
 *  Intersect two communicators                                          *
 ************************************************************************/
#ifdef USE_MPI
static inline void MPI_Group_free2(MPI_Group *group) {
    if (*group != MPI_GROUP_EMPTY) {
        // MPICH is fine with free'ing an empty group, OpenMPI crashes
        MPI_Group_free(group);
    }
}
MPI_CLASS MPI_CLASS::intersect(const MPI_CLASS &comm1, const MPI_CLASS &comm2) {
    MPI_Group group1 = MPI_GROUP_EMPTY, group2 = MPI_GROUP_EMPTY;
    if (!comm1.isNull()) {
        MPI_Group_free2(&group1);
        MPI_Comm_group(comm1.communicator, &group1);
    }
    if (!comm2.isNull()) {
        MPI_Group_free2(&group2);
        MPI_Comm_group(comm2.communicator, &group2);
    }
    MPI_Group group12;
    MPI_Group_intersection(group1, group2, &group12);
    int compare1, compare2;
    MPI_Group_compare(group1, group12, &compare1);
    MPI_Group_compare(group2, group12, &compare2);
    MPI_CLASS new_comm(MPI_CLASS_COMM_NULL);
    int size;
    MPI_Group_size(group12, &size);
    if (compare1 != MPI_UNEQUAL && size != 0) {
        // The intersection matches comm1
        new_comm = comm1;
    } else if (compare2 != MPI_UNEQUAL && size != 0) {
        // The intersection matches comm2
        new_comm = comm2;
    } else if (comm1.isNull()) {
        // comm1 is null, we can return safely (comm1 is needed for communication)
    } else {
        // The intersection is smaller than comm1 or comm2
        // Check if the new comm is nullptr for all processors
        int max_size = 0;
        MPI_Allreduce(&size, &max_size, 1, MPI_INT, MPI_MAX,
                      comm1.communicator);
        if (max_size == 0) {
            // We are dealing with completely disjoint sets
            new_comm = MPI_CLASS(MPI_CLASS_COMM_NULL, false);
        } else {
            // Create the new comm
            // Note: OpenMPI crashes if the intersection group is EMPTY for any processors
            // We will set it to SELF for the EMPTY processors, then create a nullptr comm later
            if (group12 == MPI_GROUP_EMPTY) {
                MPI_Group_free2(&group12);
                MPI_Comm_group(MPI_COMM_SELF, &group12);
            }
            MPI_Comm new_MPI_comm;
            MPI_Comm_create(comm1.communicator, group12, &new_MPI_comm);
            if (size > 0) {
                // This is the valid case where we create a new intersection comm
                new_comm = MPI_CLASS(new_MPI_comm, true);
            } else {
                // We actually want a null comm for this communicator
                new_comm = MPI_CLASS(MPI_CLASS_COMM_NULL, false);
                MPI_Comm_free(&new_MPI_comm);
            }
        }
    }
    MPI_Group_free2(&group1);
    MPI_Group_free2(&group2);
    MPI_Group_free2(&group12);
    return new_comm;
}
#else
MPI_CLASS MPI_CLASS::intersect(const MPI_CLASS &comm1, const MPI_CLASS &comm2) {
    if (comm1.isNull() || comm2.isNull())
        return MPI_CLASS(MPI_CLASS_COMM_NULL, false);
    MPI_ASSERT(comm1.comm_size == 1 && comm2.comm_size == 1);
    return comm1;
}
#endif

/************************************************************************
 *  Split a comm						                                    *
 ************************************************************************/
MPI_CLASS MPI_CLASS::split(int color, int key) const {
    if (d_isNull) {
        return MPI_CLASS(MPI_CLASS_COMM_NULL);
    } else if (comm_size == 1) {
        if (color == -1)
            return MPI_CLASS(MPI_CLASS_COMM_NULL);
        return dup();
    }
    MPI_Comm new_MPI_comm = MPI_CLASS_COMM_NULL;
#ifdef USE_MPI
    // USE MPI to split the communicator
    if (color == -1) {
        check_MPI(
            MPI_Comm_split(communicator, MPI_UNDEFINED, key, &new_MPI_comm));
    } else {
        check_MPI(MPI_Comm_split(communicator, color, key, &new_MPI_comm));
    }
#endif
    // Create the new object
    NULL_USE(key);
    MPI_CLASS new_comm(new_MPI_comm, true);
    new_comm.d_call_abort = d_call_abort;
    return new_comm;
}
MPI_CLASS MPI_CLASS::splitByNode(int key) const {
    // Check if we are dealing with a single processor (trivial case)
    if (comm_size == 1)
        return this->split(0, 0);
    // Get the node name
    std::string name = MPI_CLASS::getNodeName();
    // Gather the names from all ranks
    std::vector<std::string> list(comm_size);
    allGather(name, &list[0]);
    // Create the colors
    std::vector<int> color(comm_size, -1);
    color[0] = 0;
    for (int i = 1; i < comm_size; i++) {
        const std::string tmp1 = list[i];
        for (int j = 0; j < i; j++) {
            const std::string tmp2 = list[j];
            if (tmp1 == tmp2) {
                color[i] = color[j];
                break;
            }
            color[i] = color[i - 1] + 1;
        }
    }
    MPI_CLASS new_comm = this->split(color[comm_rank], key);
    return new_comm;
}

/************************************************************************
 *  Duplicate an exisiting comm object                                   *
 ************************************************************************/
MPI_CLASS MPI_CLASS::dup() const {
    if (d_isNull)
        return MPI_CLASS(MPI_CLASS_COMM_NULL);
    MPI_Comm new_MPI_comm = communicator;
#if defined(USE_MPI) || defined(USE_PETSC)
    // USE MPI to duplicate the communicator
    MPI_Comm_dup(communicator, &new_MPI_comm);
#else
    new_MPI_comm = uniqueGlobalComm;
    uniqueGlobalComm++;
#endif
    // Create the new comm object
    MPI_CLASS new_comm(new_MPI_comm, true);
    new_comm.d_isNull = d_isNull;
    new_comm.d_call_abort = d_call_abort;
    return new_comm;
}

/************************************************************************
 *  Get the node name                                                    *
 ************************************************************************/
std::string MPI_CLASS::getNodeName() {
#ifdef USE_MPI
    int length;
    char name[MPI_MAX_PROCESSOR_NAME + 1];
    memset(name, 0, MPI_MAX_PROCESSOR_NAME + 1);
    MPI_Get_processor_name(name, &length);
    return std::string(name);
#else
    return "Node0";
#endif
}

/************************************************************************
 *  Overload operator ==                                                 *
 ************************************************************************/
bool MPI_CLASS::operator==(const MPI_CLASS &comm) const {
    return communicator == comm.communicator;
}

/************************************************************************
 *  Overload operator !=                                                 *
 ************************************************************************/
bool MPI_CLASS::operator!=(const MPI_CLASS &comm) const {
    return communicator != comm.communicator;
}

/************************************************************************
 *  Overload operator <                                                  *
 ************************************************************************/
bool MPI_CLASS::operator<(const MPI_CLASS &comm) const {
    MPI_ASSERT(!this->d_isNull && !comm.d_isNull);
    bool flag = true;
    // First check if either communicator is NULL
    if (this->d_isNull)
        return false;
    if (comm.d_isNull)
        flag = false;
    // Use compare to check if the comms are equal
    if (compare(comm) != 0)
        return false;
    // Check that the size of the other communicator is > the current communicator size
    if (comm_size >= comm.comm_size)
        flag = false;
// Check the union of the communicator groups
// this is < comm iff this group is a subgroup of comm's group
#ifdef USE_MPI
    MPI_Group group1 = MPI_GROUP_EMPTY, group2 = MPI_GROUP_EMPTY,
              group12 = MPI_GROUP_EMPTY;
    if (!d_isNull)
        MPI_Comm_group(communicator, &group1);
    if (!comm.d_isNull)
        MPI_Comm_group(comm.communicator, &group2);
    MPI_Group_union(group1, group2, &group12);
    int compare;
    MPI_Group_compare(group2, group12, &compare);
    if (compare == MPI_UNEQUAL)
        flag = false;
    MPI_Group_free(&group1);
    MPI_Group_free(&group2);
    MPI_Group_free(&group12);
#endif
    // Perform a global reduce of the flag (equivalent to all operation)
    return allReduce(flag);
}

/************************************************************************
 *  Overload operator <=                                                 *
 ************************************************************************/
bool MPI_CLASS::operator<=(const MPI_CLASS &comm) const {
    MPI_ASSERT(!this->d_isNull && !comm.d_isNull);
    bool flag = true;
    // First check if either communicator is NULL
    if (this->d_isNull)
        return false;
    if (comm.d_isNull)
        flag = false;
#ifdef USE_MPI
    int world_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if (comm.getSize() == world_size)
        return true;
    if (getSize() == 1 && !comm.d_isNull)
        return true;
#endif
    // Use compare to check if the comms are equal
    if (compare(comm) != 0)
        return true;
    // Check that the size of the other communicator is > the current communicator size
    // this is <= comm iff this group is a subgroup of comm's group
    if (comm_size > comm.comm_size)
        flag = false;
// Check the unnion of the communicator groups
#ifdef USE_MPI
    MPI_Group group1, group2, group12;
    MPI_Comm_group(communicator, &group1);
    MPI_Comm_group(comm.communicator, &group2);
    MPI_Group_union(group1, group2, &group12);
    int compare;
    MPI_Group_compare(group2, group12, &compare);
    if (compare == MPI_UNEQUAL)
        flag = false;
    MPI_Group_free(&group1);
    MPI_Group_free(&group2);
    MPI_Group_free(&group12);
#endif
    // Perform a global reduce of the flag (equivalent to all operation)
    return allReduce(flag);
}

/************************************************************************
 *  Overload operator >                                                  *
 ************************************************************************/
bool MPI_CLASS::operator>(const MPI_CLASS &comm) const {
    bool flag = true;
    // First check if either communicator is NULL
    if (this->d_isNull)
        return false;
    if (comm.d_isNull)
        flag = false;
    // Use compare to check if the comms are equal
    if (compare(comm) != 0)
        return false;
    // Check that the size of the other communicator is > the current communicator size
    if (comm_size <= comm.comm_size)
        flag = false;
// Check the unnion of the communicator groups
// this is > comm iff comm's group is a subgroup of this group
#ifdef USE_MPI
    MPI_Group group1 = MPI_GROUP_EMPTY, group2 = MPI_GROUP_EMPTY,
              group12 = MPI_GROUP_EMPTY;
    if (!d_isNull)
        MPI_Comm_group(communicator, &group1);
    if (!comm.d_isNull)
        MPI_Comm_group(comm.communicator, &group2);
    MPI_Group_union(group1, group2, &group12);
    int compare;
    MPI_Group_compare(group1, group12, &compare);
    if (compare == MPI_UNEQUAL)
        flag = false;
    MPI_Group_free(&group1);
    MPI_Group_free(&group2);
    MPI_Group_free(&group12);
#endif
    // Perform a global reduce of the flag (equivalent to all operation)
    return allReduce(flag);
}

/************************************************************************
 *  Overload operator >=                                                 *
 ************************************************************************/
bool MPI_CLASS::operator>=(const MPI_CLASS &comm) const {
    bool flag = true;
    // First check if either communicator is NULL
    if (this->d_isNull)
        return false;
    if (comm.d_isNull)
        flag = false;
#ifdef USE_MPI
    int world_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if (getSize() == world_size)
        return true;
    if (comm.getSize() == 1 && !comm.d_isNull)
        return true;
#endif
    // Use compare to check if the comms are equal
    if (compare(comm) != 0)
        return true;
    // Check that the size of the other communicator is > the current communicator size
    if (comm_size < comm.comm_size)
        flag = false;
// Check the unnion of the communicator groups
// this is >= comm iff comm's group is a subgroup of this group
#ifdef USE_MPI
    MPI_Group group1 = MPI_GROUP_EMPTY, group2 = MPI_GROUP_EMPTY,
              group12 = MPI_GROUP_EMPTY;
    if (!d_isNull)
        MPI_Comm_group(communicator, &group1);
    if (!comm.d_isNull)
        MPI_Comm_group(comm.communicator, &group2);
    MPI_Group_union(group1, group2, &group12);
    int compare;
    MPI_Group_compare(group1, group12, &compare);
    if (compare == MPI_UNEQUAL)
        flag = false;
    MPI_Group_free(&group1);
    MPI_Group_free(&group2);
    MPI_Group_free(&group12);
#endif
    // Perform a global reduce of the flag (equivalent to all operation)
    return allReduce(flag);
}

/************************************************************************
 *  Compare two comm objects                                             *
 ************************************************************************/
int MPI_CLASS::compare(const MPI_CLASS &comm) const {
    if (communicator == comm.communicator)
        return 1;
#ifdef USE_MPI
    if (d_isNull || comm.d_isNull)
        return 0;
    int result;
    check_MPI(MPI_Comm_compare(communicator, comm.communicator, &result));
    if (result == MPI_IDENT)
        return 2;
    else if (result == MPI_CONGRUENT)
        return 3;
    else if (result == MPI_SIMILAR)
        return 4;
    else if (result == MPI_UNEQUAL)
        return 0;
    MPI_ERROR("Unknown results from comm compare");
#else
    if (comm.communicator == MPI_COMM_NULL || communicator == MPI_COMM_NULL)
        return 0;
    else
        return 3;
#endif
    return 0;
}

/************************************************************************
 *  Abort the program.                                                   *
 ************************************************************************/
void MPI_CLASS::setCallAbortInSerialInsteadOfExit(bool flag) {
    d_call_abort = flag;
}
void MPI_CLASS::abort() const {
#ifdef USE_MPI
    MPI_Comm comm = communicator;
    if (comm == MPI_COMM_NULL)
        comm = MPI_COMM_WORLD;
    if (!MPI_active()) {
        // MPI is not availible
        exit(-1);
    } else if (comm_size > 1) {
        MPI_Abort(comm, -1);
    } else if (d_call_abort) {
        MPI_Abort(comm, -1);
    } else {
        exit(-1);
    }
#else
    exit(-1);
#endif
}

/************************************************************************
 *  newTag                                                               *
 ************************************************************************/
int MPI_CLASS::newTag() {
#ifdef USE_MPI
    // Syncronize the processes to ensure all ranks enter this call
    // Needed so the count will match
    barrier();
    // Return and increment the tag
    int tag = (*d_currentTag)++;
    MPI_INSIST(tag <= d_maxTag, "Maximum number of tags exceeded\n");
    return tag;
#else
    static int globalCurrentTag = 1;
    return globalCurrentTag++;
#endif
}

/************************************************************************
 *  allReduce                                                            *
 ************************************************************************/
bool MPI_CLASS::allReduce(const bool value) const {
    bool ret = value;
    if (comm_size > 1) {
#ifdef USE_MPI
        MPI_Allreduce((void *)&value, (void *)&ret, 1, MPI_UNSIGNED_CHAR,
                      MPI_MIN, communicator);
#else
        MPI_ERROR("This shouldn't be possible");
#endif
    }
    return ret;
}

/************************************************************************
 *  anyReduce                                                            *
 ************************************************************************/
bool MPI_CLASS::anyReduce(const bool value) const {
    bool ret = value;
    if (comm_size > 1) {
#ifdef USE_MPI
        MPI_Allreduce((void *)&value, (void *)&ret, 1, MPI_UNSIGNED_CHAR,
                      MPI_MAX, communicator);
#else
        MPI_ERROR("This shouldn't be possible");
#endif
    }
    return ret;
}

/************************************************************************
 *  call_sumReduce                                                       *
 *  Note: these specializations are only called when using MPI.          *
 ************************************************************************/
#ifdef USE_MPI
// unsigned char
template <>
void MPI_CLASS::call_sumReduce<unsigned char>(const unsigned char *send,
                                              unsigned char *recv,
                                              int n) const {
    PROFILE_START("sumReduce1<unsigned char>", profile_level);
    MPI_Allreduce((void *)send, (void *)recv, n, MPI_UNSIGNED_CHAR, MPI_SUM,
                  communicator);
    PROFILE_STOP("sumReduce1<unsigned char>", profile_level);
}
template <>
void MPI_CLASS::call_sumReduce<unsigned char>(unsigned char *x, int n) const {
    PROFILE_START("sumReduce2<unsigned char>", profile_level);
    auto send = x;
    auto recv = new unsigned char[n];
    MPI_Allreduce(send, recv, n, MPI_UNSIGNED_CHAR, MPI_SUM, communicator);
    for (int i = 0; i < n; i++)
        x[i] = recv[i];
    delete[] recv;
    PROFILE_STOP("sumReduce2<unsigned char>", profile_level);
}
// char
template <>
void MPI_CLASS::call_sumReduce<char>(const char *send, char *recv,
                                     int n) const {
    PROFILE_START("sumReduce1<char>", profile_level);
    MPI_Allreduce((void *)send, (void *)recv, n, MPI_SIGNED_CHAR, MPI_SUM,
                  communicator);
    PROFILE_STOP("sumReduce1<char>", profile_level);
}
template <> void MPI_CLASS::call_sumReduce<char>(char *x, int n) const {
    PROFILE_START("sumReduce2<char>", profile_level);
    auto send = x;
    auto recv = new char[n];
    MPI_Allreduce(send, recv, n, MPI_SIGNED_CHAR, MPI_SUM, communicator);
    for (int i = 0; i < n; i++)
        x[i] = recv[i];
    delete[] recv;
    PROFILE_STOP("sumReduce2<char>", profile_level);
}
// unsigned int
template <>
void MPI_CLASS::call_sumReduce<unsigned int>(const unsigned int *send,
                                             unsigned int *recv, int n) const {
    PROFILE_START("sumReduce1<unsigned int>", profile_level);
    MPI_Allreduce((void *)send, (void *)recv, n, MPI_UNSIGNED, MPI_SUM,
                  communicator);
    PROFILE_STOP("sumReduce1<unsigned int>", profile_level);
}
template <>
void MPI_CLASS::call_sumReduce<unsigned int>(unsigned int *x, int n) const {
    PROFILE_START("sumReduce2<unsigned int>", profile_level);
    auto send = x;
    auto recv = new unsigned int[n];
    MPI_Allreduce(send, recv, n, MPI_UNSIGNED, MPI_SUM, communicator);
    for (int i = 0; i < n; i++)
        x[i] = recv[i];
    delete[] recv;
    PROFILE_STOP("sumReduce2<unsigned int>", profile_level);
}
// int
template <>
void MPI_CLASS::call_sumReduce<int>(const int *send, int *recv, int n) const {
    PROFILE_START("sumReduce1<int>", profile_level);
    MPI_Allreduce((void *)send, (void *)recv, n, MPI_INT, MPI_SUM,
                  communicator);
    PROFILE_STOP("sumReduce1<int>", profile_level);
}
template <> void MPI_CLASS::call_sumReduce<int>(int *x, int n) const {
    PROFILE_START("sumReduce2<int>", profile_level);
    auto send = x;
    auto recv = new int[n];
    MPI_Allreduce(send, recv, n, MPI_INT, MPI_SUM, communicator);
    for (int i = 0; i < n; i++)
        x[i] = recv[i];
    delete[] recv;
    PROFILE_STOP("sumReduce2<int>", profile_level);
}
// long int
template <>
void MPI_CLASS::call_sumReduce<long int>(const long int *send, long int *recv,
                                         int n) const {
    PROFILE_START("sumReduce1<long int>", profile_level);
    MPI_Allreduce((void *)send, (void *)recv, n, MPI_LONG, MPI_SUM,
                  communicator);
    PROFILE_STOP("sumReduce1<long int>", profile_level);
}
template <> void MPI_CLASS::call_sumReduce<long int>(long int *x, int n) const {
    PROFILE_START("sumReduce2<long int>", profile_level);
    auto send = x;
    auto recv = new long int[n];
    MPI_Allreduce(send, recv, n, MPI_LONG, MPI_SUM, communicator);
    for (int i = 0; i < n; i++)
        x[i] = recv[i];
    delete[] recv;
    PROFILE_STOP("sumReduce2<long int>", profile_level);
}
// unsigned long int
template <>
void MPI_CLASS::call_sumReduce<unsigned long>(const unsigned long *send,
                                              unsigned long *recv,
                                              int n) const {
    PROFILE_START("sumReduce1<unsigned long>", profile_level);
    MPI_Allreduce((void *)send, (void *)recv, n, MPI_UNSIGNED_LONG, MPI_SUM,
                  communicator);
    PROFILE_STOP("sumReduce1<unsigned long>", profile_level);
}
template <>
void MPI_CLASS::call_sumReduce<unsigned long>(unsigned long *x, int n) const {
    PROFILE_START("sumReduce2<unsigned long>", profile_level);
    auto send = x;
    auto recv = new unsigned long int[n];
    MPI_Allreduce(send, recv, n, MPI_UNSIGNED_LONG, MPI_SUM, communicator);
    for (int i = 0; i < n; i++)
        x[i] = recv[i];
    delete[] recv;
    PROFILE_STOP("sumReduce2<unsigned long>", profile_level);
}
// size_t
#ifdef USE_WINDOWS
template <>
void MPI_CLASS::call_sumReduce<size_t>(const size_t *send, size_t *recv,
                                       int n) const {
    MPI_ASSERT(MPI_SIZE_T != 0);
    PROFILE_START("sumReduce1<size_t>", profile_level);
    MPI_Allreduce((void *)send, (void *)recv, n, MPI_SIZE_T, MPI_SUM,
                  communicator);
    PROFILE_STOP("sumReduce1<size_t>", profile_level);
}
template <> void MPI_CLASS::call_sumReduce<size_t>(size_t *x, int n) const {
    MPI_ASSERT(MPI_SIZE_T != 0);
    PROFILE_START("sumReduce2<size_t>", profile_level);
    auto send = x;
    auto recv = new size_t[n];
    MPI_Allreduce((void *)send, (void *)recv, n, MPI_SIZE_T, MPI_SUM,
                  communicator);
    for (int i = 0; i < n; i++)
        x[i] = recv[i];
    delete[] recv;
    PROFILE_STOP("sumReduce2<size_t>", profile_level);
}
#endif
// float
template <>
void MPI_CLASS::call_sumReduce<float>(const float *send, float *recv,
                                      int n) const {
    PROFILE_START("sumReduce1<float>", profile_level);
    MPI_Allreduce((void *)send, (void *)recv, n, MPI_FLOAT, MPI_SUM,
                  communicator);
    PROFILE_STOP("sumReduce1<float>", profile_level);
}
template <> void MPI_CLASS::call_sumReduce<float>(float *x, int n) const {
    PROFILE_START("sumReduce2<float>", profile_level);
    auto send = x;
    auto recv = new float[n];
    MPI_Allreduce(send, recv, n, MPI_FLOAT, MPI_SUM, communicator);
    for (int i = 0; i < n; i++)
        x[i] = recv[i];
    delete[] recv;
    PROFILE_STOP("sumReduce2<float>", profile_level);
}
// double
template <>
void MPI_CLASS::call_sumReduce<double>(const double *send, double *recv,
                                       int n) const {
    PROFILE_START("sumReduce1<double>", profile_level);
    MPI_Allreduce((void *)send, (void *)recv, n, MPI_DOUBLE, MPI_SUM,
                  communicator);
    PROFILE_STOP("sumReduce1<double>", profile_level);
}
template <> void MPI_CLASS::call_sumReduce<double>(double *x, int n) const {
    PROFILE_START("sumReduce2<double>", profile_level);
    auto send = x;
    auto recv = new double[n];
    MPI_Allreduce(send, recv, n, MPI_DOUBLE, MPI_SUM, communicator);
    for (int i = 0; i < n; i++)
        x[i] = recv[i];
    delete[] recv;
    PROFILE_STOP("sumReduce2<double>", profile_level);
}
// std::complex<double>
template <>
void MPI_CLASS::call_sumReduce<std::complex<double>>(
    const std::complex<double> *x, std::complex<double> *y, int n) const {
    PROFILE_START("sumReduce1<complex double>", profile_level);
    auto send = new double[2 * n];
    auto recv = new double[2 * n];
    for (int i = 0; i < n; i++) {
        send[2 * i + 0] = real(x[i]);
        send[2 * i + 1] = imag(x[i]);
    }
    MPI_Allreduce((void *)send, (void *)recv, 2 * n, MPI_DOUBLE, MPI_SUM,
                  communicator);
    for (int i = 0; i < n; i++)
        y[i] = std::complex<double>(recv[2 * i + 0], recv[2 * i + 1]);
    delete[] send;
    delete[] recv;
    PROFILE_STOP("sumReduce1<complex double>", profile_level);
}
template <>
void MPI_CLASS::call_sumReduce<std::complex<double>>(std::complex<double> *x,
                                                     int n) const {
    PROFILE_START("sumReduce2<complex double>", profile_level);
    auto send = new double[2 * n];
    auto recv = new double[2 * n];
    for (int i = 0; i < n; i++) {
        send[2 * i + 0] = real(x[i]);
        send[2 * i + 1] = imag(x[i]);
    }
    MPI_Allreduce(send, recv, 2 * n, MPI_DOUBLE, MPI_SUM, communicator);
    for (int i = 0; i < n; i++)
        x[i] = std::complex<double>(recv[2 * i + 0], recv[2 * i + 1]);
    delete[] send;
    delete[] recv;
    PROFILE_STOP("sumReduce2<complex double>", profile_level);
}
#endif

/************************************************************************
 *  call_minReduce                                                       *
 *  Note: these specializations are only called when using MPI.          *
 ************************************************************************/
#ifdef USE_MPI
// unsigned char
template <>
void MPI_CLASS::call_minReduce<unsigned char>(const unsigned char *send,
                                              unsigned char *recv, int n,
                                              int *comm_rank_of_min) const {
    if (comm_rank_of_min == nullptr) {
        PROFILE_START("minReduce1<unsigned char>", profile_level);
        MPI_Allreduce((void *)send, (void *)recv, n, MPI_UNSIGNED_CHAR, MPI_MIN,
                      communicator);
        PROFILE_STOP("minReduce1<unsigned char>", profile_level);
    } else {
        auto tmp = new int[n];
        for (int i = 0; i < n; i++)
            tmp[i] = send[i];
        call_minReduce<int>(tmp, n, comm_rank_of_min);
        for (int i = 0; i < n; i++)
            recv[i] = static_cast<unsigned char>(tmp[i]);
        delete[] tmp;
    }
}
template <>
void MPI_CLASS::call_minReduce<unsigned char>(unsigned char *x, int n,
                                              int *comm_rank_of_min) const {
    if (comm_rank_of_min == nullptr) {
        PROFILE_START("minReduce2<unsigned char>", profile_level);
        auto send = x;
        auto recv = new unsigned char[n];
        MPI_Allreduce(send, recv, n, MPI_UNSIGNED_CHAR, MPI_MIN, communicator);
        for (int i = 0; i < n; i++)
            x[i] = recv[i];
        delete[] recv;
        PROFILE_STOP("minReduce2<unsigned char>", profile_level);
    } else {
        auto tmp = new int[n];
        for (int i = 0; i < n; i++)
            tmp[i] = x[i];
        call_minReduce<int>(tmp, n, comm_rank_of_min);
        for (int i = 0; i < n; i++)
            x[i] = static_cast<unsigned char>(tmp[i]);
        delete[] tmp;
    }
}
// char
template <>
void MPI_CLASS::call_minReduce<char>(const char *send, char *recv, int n,
                                     int *comm_rank_of_min) const {
    if (comm_rank_of_min == nullptr) {
        PROFILE_START("minReduce1<char>", profile_level);
        MPI_Allreduce((void *)send, (void *)recv, n, MPI_SIGNED_CHAR, MPI_MIN,
                      communicator);
        PROFILE_STOP("minReduce1<char>", profile_level);
    } else {
        auto tmp = new int[n];
        for (int i = 0; i < n; i++)
            tmp[i] = send[i];
        call_minReduce<int>(tmp, n, comm_rank_of_min);
        for (int i = 0; i < n; i++)
            recv[i] = static_cast<char>(tmp[i]);
        delete[] tmp;
    }
}
template <>
void MPI_CLASS::call_minReduce<char>(char *x, int n,
                                     int *comm_rank_of_min) const {
    if (comm_rank_of_min == nullptr) {
        PROFILE_START("minReduce2<char>", profile_level);
        auto send = x;
        auto recv = new char[n];
        MPI_Allreduce(send, recv, n, MPI_SIGNED_CHAR, MPI_MIN, communicator);
        for (int i = 0; i < n; i++)
            x[i] = recv[i];
        delete[] recv;
        PROFILE_STOP("minReduce2<char>", profile_level);
    } else {
        auto tmp = new int[n];
        for (int i = 0; i < n; i++)
            tmp[i] = x[i];
        call_minReduce<int>(tmp, n, comm_rank_of_min);
        for (int i = 0; i < n; i++)
            x[i] = static_cast<char>(tmp[i]);
        delete[] tmp;
    }
}
// unsigned int
template <>
void MPI_CLASS::call_minReduce<unsigned int>(const unsigned int *send,
                                             unsigned int *recv, int n,
                                             int *comm_rank_of_min) const {
    if (comm_rank_of_min == nullptr) {
        PROFILE_START("minReduce1<unsigned int>", profile_level);
        MPI_Allreduce((void *)send, (void *)recv, n, MPI_UNSIGNED, MPI_MIN,
                      communicator);
        PROFILE_STOP("minReduce1<unsigned int>", profile_level);
    } else {
        auto tmp = new int[n];
        for (int i = 0; i < n; i++)
            tmp[i] = unsigned_to_signed(send[i]);
        call_minReduce<int>(tmp, n, comm_rank_of_min);
        for (int i = 0; i < n; i++)
            recv[i] = signed_to_unsigned(tmp[i]);
        delete[] tmp;
    }
}
template <>
void MPI_CLASS::call_minReduce<unsigned int>(unsigned int *x, int n,
                                             int *comm_rank_of_min) const {
    if (comm_rank_of_min == nullptr) {
        PROFILE_START("minReduce2<unsigned int>", profile_level);
        auto send = x;
        auto recv = new unsigned int[n];
        MPI_Allreduce(send, recv, n, MPI_UNSIGNED, MPI_MIN, communicator);
        for (int i = 0; i < n; i++)
            x[i] = recv[i];
        delete[] recv;
        PROFILE_STOP("minReduce2<unsigned int>", profile_level);
    } else {
        auto tmp = new int[n];
        for (int i = 0; i < n; i++)
            tmp[i] = unsigned_to_signed(x[i]);
        call_minReduce<int>(tmp, n, comm_rank_of_min);
        for (int i = 0; i < n; i++)
            x[i] = signed_to_unsigned(tmp[i]);
        delete[] tmp;
    }
}
// int
template <>
void MPI_CLASS::call_minReduce<int>(const int *x, int *y, int n,
                                    int *comm_rank_of_min) const {
    PROFILE_START("minReduce1<int>", profile_level);
    if (comm_rank_of_min == nullptr) {
        MPI_Allreduce((void *)x, (void *)y, n, MPI_INT, MPI_MIN, communicator);
    } else {
        auto recv = new IntIntStruct[n];
        auto send = new IntIntStruct[n];
        for (int i = 0; i < n; ++i) {
            send[i].j = x[i];
            send[i].i = comm_rank;
        }
        MPI_Allreduce(send, recv, n, MPI_2INT, MPI_MINLOC, communicator);
        for (int i = 0; i < n; ++i) {
            y[i] = recv[i].j;
            comm_rank_of_min[i] = recv[i].i;
        }
        delete[] recv;
        delete[] send;
    }
    PROFILE_STOP("minReduce1<int>", profile_level);
}
template <>
void MPI_CLASS::call_minReduce<int>(int *x, int n,
                                    int *comm_rank_of_min) const {
    PROFILE_START("minReduce2<int>", profile_level);
    if (comm_rank_of_min == nullptr) {
        auto send = x;
        auto recv = new int[n];
        MPI_Allreduce(send, recv, n, MPI_INT, MPI_MIN, communicator);
        for (int i = 0; i < n; i++)
            x[i] = recv[i];
        delete[] recv;
    } else {
        auto recv = new IntIntStruct[n];
        auto send = new IntIntStruct[n];
        for (int i = 0; i < n; ++i) {
            send[i].j = x[i];
            send[i].i = comm_rank;
        }
        MPI_Allreduce(send, recv, n, MPI_2INT, MPI_MINLOC, communicator);
        for (int i = 0; i < n; ++i) {
            x[i] = recv[i].j;
            comm_rank_of_min[i] = recv[i].i;
        }
        delete[] recv;
        delete[] send;
    }
    PROFILE_STOP("minReduce2<int>", profile_level);
}
// unsigned long int
template <>
void MPI_CLASS::call_minReduce<unsigned long int>(const unsigned long int *send,
                                                  unsigned long int *recv,
                                                  int n,
                                                  int *comm_rank_of_min) const {
    if (comm_rank_of_min == nullptr) {
        PROFILE_START("minReduce1<unsigned long>", profile_level);
        MPI_Allreduce((void *)send, (void *)recv, n, MPI_UNSIGNED_LONG, MPI_MIN,
                      communicator);
        PROFILE_STOP("minReduce1<unsigned long>", profile_level);
    } else {
        auto tmp = new long int[n];
        for (int i = 0; i < n; i++)
            tmp[i] = unsigned_to_signed(send[i]);
        call_minReduce<long int>(tmp, n, comm_rank_of_min);
        for (int i = 0; i < n; i++)
            recv[i] = signed_to_unsigned(tmp[i]);
        delete[] tmp;
    }
}
template <>
void MPI_CLASS::call_minReduce<unsigned long int>(unsigned long int *x, int n,
                                                  int *comm_rank_of_min) const {
    if (comm_rank_of_min == nullptr) {
        PROFILE_START("minReduce2<unsigned long>", profile_level);
        auto send = x;
        auto recv = new unsigned long int[n];
        MPI_Allreduce(send, recv, n, MPI_UNSIGNED_LONG, MPI_MIN, communicator);
        for (int i = 0; i < n; i++)
            x[i] = recv[i];
        delete[] recv;
        PROFILE_STOP("minReduce2<unsigned long>", profile_level);
    } else {
        auto tmp = new long int[n];
        for (int i = 0; i < n; i++)
            tmp[i] = unsigned_to_signed(x[i]);
        call_minReduce<long int>(tmp, n, comm_rank_of_min);
        for (int i = 0; i < n; i++)
            x[i] = signed_to_unsigned(tmp[i]);
        delete[] tmp;
    }
}
// long int
template <>
void MPI_CLASS::call_minReduce<long int>(const long int *x, long int *y, int n,
                                         int *comm_rank_of_min) const {
    PROFILE_START("minReduce1<long int>", profile_level);
    if (comm_rank_of_min == nullptr) {
        MPI_Allreduce((void *)x, (void *)y, n, MPI_LONG, MPI_MIN, communicator);
    } else {
        auto recv = new LongIntStruct[n];
        auto send = new LongIntStruct[n];
        for (int i = 0; i < n; ++i) {
            send[i].j = x[i];
            send[i].i = comm_rank;
        }
        MPI_Allreduce(send, recv, n, MPI_LONG_INT, MPI_MINLOC, communicator);
        for (int i = 0; i < n; ++i) {
            y[i] = recv[i].j;
            comm_rank_of_min[i] = recv[i].i;
        }
        delete[] recv;
        delete[] send;
    }
    PROFILE_STOP("minReduce1<long int>", profile_level);
}
template <>
void MPI_CLASS::call_minReduce<long int>(long int *x, int n,
                                         int *comm_rank_of_min) const {
    PROFILE_START("minReduce2<long int>", profile_level);
    if (comm_rank_of_min == nullptr) {
        auto send = x;
        auto recv = new long int[n];
        MPI_Allreduce(send, recv, n, MPI_LONG, MPI_MIN, communicator);
        for (long int i = 0; i < n; i++)
            x[i] = recv[i];
        delete[] recv;
    } else {
        auto recv = new LongIntStruct[n];
        auto send = new LongIntStruct[n];
        for (int i = 0; i < n; ++i) {
            send[i].j = x[i];
            send[i].i = comm_rank;
        }
        MPI_Allreduce(send, recv, n, MPI_LONG_INT, MPI_MINLOC, communicator);
        for (int i = 0; i < n; ++i) {
            x[i] = recv[i].j;
            comm_rank_of_min[i] = recv[i].i;
        }
        delete[] recv;
        delete[] send;
    }
    PROFILE_STOP("minReduce2<long int>", profile_level);
}
// unsigned long long int
template <>
void MPI_CLASS::call_minReduce<unsigned long long int>(
    const unsigned long long int *send, unsigned long long int *recv, int n,
    int *comm_rank_of_min) const {
    PROFILE_START("minReduce1<long int>", profile_level);
    if (comm_rank_of_min == nullptr) {
        auto x = new long long int[n];
        auto y = new long long int[n];
        for (int i = 0; i < n; i++)
            x[i] = unsigned_to_signed(send[i]);
        MPI_Allreduce((void *)x, (void *)y, n, MPI_LONG_LONG_INT, MPI_MIN,
                      communicator);
        for (int i = 0; i < n; i++)
            recv[i] = signed_to_unsigned(y[i]);
        delete[] x;
        delete[] y;
    } else {
        printf("minReduce<long long int> will use double\n");
        auto tmp = new double[n];
        for (int i = 0; i < n; i++)
            tmp[i] = static_cast<double>(send[i]);
        call_minReduce<double>(tmp, n, comm_rank_of_min);
        for (int i = 0; i < n; i++)
            recv[i] = static_cast<long long int>(tmp[i]);
        delete[] tmp;
    }
    PROFILE_STOP("minReduce1<long int>", profile_level);
}
template <>
void MPI_CLASS::call_minReduce<unsigned long long int>(
    unsigned long long int *x, int n, int *comm_rank_of_min) const {
    auto recv = new unsigned long long int[n];
    call_minReduce<unsigned long long int>(x, recv, n, comm_rank_of_min);
    for (int i = 0; i < n; i++)
        x[i] = recv[i];
    delete[] recv;
}
// long long int
template <>
void MPI_CLASS::call_minReduce<long long int>(const long long int *x,
                                              long long int *y, int n,
                                              int *comm_rank_of_min) const {
    PROFILE_START("minReduce1<long int>", profile_level);
    if (comm_rank_of_min == nullptr) {
        MPI_Allreduce((void *)x, (void *)y, n, MPI_LONG_LONG_INT, MPI_MIN,
                      communicator);
    } else {
        printf("minReduce<long long int> will use double\n");
        auto tmp = new double[n];
        for (int i = 0; i < n; i++)
            tmp[i] = static_cast<double>(x[i]);
        call_minReduce<double>(tmp, n, comm_rank_of_min);
        for (int i = 0; i < n; i++)
            y[i] = static_cast<long long int>(tmp[i]);
        delete[] tmp;
    }
    PROFILE_STOP("minReduce1<long int>", profile_level);
}
template <>
void MPI_CLASS::call_minReduce<long long int>(long long int *x, int n,
                                              int *comm_rank_of_min) const {
    auto recv = new long long int[n];
    call_minReduce<long long int>(x, recv, n, comm_rank_of_min);
    for (int i = 0; i < n; i++)
        x[i] = signed_to_unsigned(recv[i]);
    delete[] recv;
}
// float
template <>
void MPI_CLASS::call_minReduce<float>(const float *x, float *y, int n,
                                      int *comm_rank_of_min) const {
    PROFILE_START("minReduce1<float>", profile_level);
    if (comm_rank_of_min == nullptr) {
        MPI_Allreduce((void *)x, (void *)y, n, MPI_INT, MPI_MIN, communicator);
    } else {
        auto recv = new FloatIntStruct[n];
        auto send = new FloatIntStruct[n];
        for (int i = 0; i < n; ++i) {
            send[i].f = x[i];
            send[i].i = comm_rank;
        }
        MPI_Allreduce(send, recv, n, MPI_FLOAT_INT, MPI_MINLOC, communicator);
        for (int i = 0; i < n; ++i) {
            y[i] = recv[i].f;
            comm_rank_of_min[i] = recv[i].i;
        }
        delete[] recv;
        delete[] send;
    }
    PROFILE_STOP("minReduce1<float>", profile_level);
}
template <>
void MPI_CLASS::call_minReduce<float>(float *x, int n,
                                      int *comm_rank_of_min) const {
    PROFILE_START("minReduce2<float>", profile_level);
    if (comm_rank_of_min == nullptr) {
        auto send = x;
        auto recv = new float[n];
        MPI_Allreduce(send, recv, n, MPI_FLOAT, MPI_MIN, communicator);
        for (int i = 0; i < n; i++)
            x[i] = recv[i];
        delete[] recv;
    } else {
        auto recv = new FloatIntStruct[n];
        auto send = new FloatIntStruct[n];
        for (int i = 0; i < n; ++i) {
            send[i].f = x[i];
            send[i].i = comm_rank;
        }
        MPI_Allreduce(send, recv, n, MPI_FLOAT_INT, MPI_MINLOC, communicator);
        for (int i = 0; i < n; ++i) {
            x[i] = recv[i].f;
            comm_rank_of_min[i] = recv[i].i;
        }
        delete[] recv;
        delete[] send;
    }
    PROFILE_STOP("minReduce2<float>", profile_level);
}
// double
template <>
void MPI_CLASS::call_minReduce<double>(const double *x, double *y, int n,
                                       int *comm_rank_of_min) const {
    PROFILE_START("minReduce1<double>", profile_level);
    if (comm_rank_of_min == nullptr) {
        MPI_Allreduce((void *)x, (void *)y, n, MPI_DOUBLE, MPI_MIN,
                      communicator);
    } else {
        auto recv = new DoubleIntStruct[n];
        auto send = new DoubleIntStruct[n];
        for (int i = 0; i < n; ++i) {
            send[i].d = x[i];
            send[i].i = comm_rank;
        }
        MPI_Allreduce(send, recv, n, MPI_DOUBLE_INT, MPI_MINLOC, communicator);
        for (int i = 0; i < n; ++i) {
            y[i] = recv[i].d;
            comm_rank_of_min[i] = recv[i].i;
        }
        delete[] recv;
        delete[] send;
    }
    PROFILE_STOP("minReduce1<double>", profile_level);
}
template <>
void MPI_CLASS::call_minReduce<double>(double *x, int n,
                                       int *comm_rank_of_min) const {
    PROFILE_START("minReduce2<double>", profile_level);
    if (comm_rank_of_min == nullptr) {
        auto send = x;
        auto recv = new double[n];
        MPI_Allreduce(send, recv, n, MPI_DOUBLE, MPI_MIN, communicator);
        for (int i = 0; i < n; i++)
            x[i] = recv[i];
        delete[] recv;
    } else {
        auto recv = new DoubleIntStruct[n];
        auto send = new DoubleIntStruct[n];
        for (int i = 0; i < n; ++i) {
            send[i].d = x[i];
            send[i].i = comm_rank;
        }
        MPI_Allreduce(send, recv, n, MPI_DOUBLE_INT, MPI_MINLOC, communicator);
        for (int i = 0; i < n; ++i) {
            x[i] = recv[i].d;
            comm_rank_of_min[i] = recv[i].i;
        }
        delete[] recv;
        delete[] send;
    }
    PROFILE_STOP("minReduce2<double>", profile_level);
}
#endif

/************************************************************************
 *  call_maxReduce                                                    *
 *  Note: these specializations are only called when using MPI.          *
 ************************************************************************/
#ifdef USE_MPI
// unsigned char
template <>
void MPI_CLASS::call_maxReduce<unsigned char>(const unsigned char *send,
                                              unsigned char *recv, int n,
                                              int *comm_rank_of_max) const {
    if (comm_rank_of_max == nullptr) {
        PROFILE_START("maxReduce1<unsigned char>", profile_level);
        MPI_Allreduce((void *)send, (void *)recv, n, MPI_UNSIGNED_CHAR, MPI_MAX,
                      communicator);
        PROFILE_STOP("maxReduce1<unsigned char>", profile_level);
    } else {
        auto tmp = new int[n];
        for (int i = 0; i < n; i++)
            tmp[i] = send[i];
        call_maxReduce<int>(tmp, n, comm_rank_of_max);
        for (int i = 0; i < n; i++)
            recv[i] = static_cast<unsigned char>(tmp[i]);
        delete[] tmp;
    }
}
template <>
void MPI_CLASS::call_maxReduce<unsigned char>(unsigned char *x, int n,
                                              int *comm_rank_of_max) const {
    if (comm_rank_of_max == nullptr) {
        PROFILE_START("maxReduce2<unsigned char>", profile_level);
        auto send = x;
        auto recv = new unsigned char[n];
        MPI_Allreduce(send, recv, n, MPI_UNSIGNED_CHAR, MPI_MAX, communicator);
        for (int i = 0; i < n; i++)
            x[i] = recv[i];
        delete[] recv;
        PROFILE_STOP("maxReduce2<unsigned char>", profile_level);
    } else {
        auto tmp = new int[n];
        for (int i = 0; i < n; i++)
            tmp[i] = x[i];
        call_maxReduce<int>(tmp, n, comm_rank_of_max);
        for (int i = 0; i < n; i++)
            x[i] = static_cast<unsigned char>(tmp[i]);
        delete[] tmp;
    }
}
// char
template <>
void MPI_CLASS::call_maxReduce<char>(const char *send, char *recv, int n,
                                     int *comm_rank_of_max) const {
    if (comm_rank_of_max == nullptr) {
        PROFILE_START("maxReduce1<char>", profile_level);
        MPI_Allreduce((void *)send, (void *)recv, n, MPI_SIGNED_CHAR, MPI_MAX,
                      communicator);
        PROFILE_STOP("maxReduce1<char>", profile_level);
    } else {
        auto tmp = new int[n];
        for (int i = 0; i < n; i++)
            tmp[i] = send[i];
        call_maxReduce<int>(tmp, n, comm_rank_of_max);
        for (int i = 0; i < n; i++)
            recv[i] = static_cast<char>(tmp[i]);
        delete[] tmp;
    }
}
template <>
void MPI_CLASS::call_maxReduce<char>(char *x, int n,
                                     int *comm_rank_of_max) const {
    if (comm_rank_of_max == nullptr) {
        PROFILE_START("maxReduce2<char>", profile_level);
        auto send = x;
        auto recv = new char[n];
        MPI_Allreduce(send, recv, n, MPI_SIGNED_CHAR, MPI_MAX, communicator);
        for (int i = 0; i < n; i++)
            x[i] = recv[i];
        delete[] recv;
        PROFILE_STOP("maxReduce2<char>", profile_level);
    } else {
        auto tmp = new int[n];
        for (int i = 0; i < n; i++)
            tmp[i] = x[i];
        call_maxReduce<int>(tmp, n, comm_rank_of_max);
        for (int i = 0; i < n; i++)
            x[i] = static_cast<char>(tmp[i]);
        delete[] tmp;
    }
}
// unsigned int
template <>
void MPI_CLASS::call_maxReduce<unsigned int>(const unsigned int *send,
                                             unsigned int *recv, int n,
                                             int *comm_rank_of_max) const {
    if (comm_rank_of_max == nullptr) {
        PROFILE_START("maxReduce1<unsigned int>", profile_level);
        MPI_Allreduce((void *)send, (void *)recv, n, MPI_UNSIGNED, MPI_MAX,
                      communicator);
        PROFILE_STOP("maxReduce1<unsigned int>", profile_level);
    } else {
        auto tmp = new int[n];
        for (int i = 0; i < n; i++)
            tmp[i] = unsigned_to_signed(send[i]);
        call_maxReduce<int>(tmp, n, comm_rank_of_max);
        for (int i = 0; i < n; i++)
            recv[i] = signed_to_unsigned(tmp[i]);
        delete[] tmp;
    }
}
template <>
void MPI_CLASS::call_maxReduce<unsigned int>(unsigned int *x, int n,
                                             int *comm_rank_of_max) const {
    if (comm_rank_of_max == nullptr) {
        PROFILE_START("maxReduce2<unsigned int>", profile_level);
        auto send = x;
        auto recv = new unsigned int[n];
        MPI_Allreduce(send, recv, n, MPI_UNSIGNED, MPI_MAX, communicator);
        for (int i = 0; i < n; i++)
            x[i] = recv[i];
        delete[] recv;
        PROFILE_STOP("maxReduce2<unsigned int>", profile_level);
    } else {
        auto tmp = new int[n];
        for (int i = 0; i < n; i++)
            tmp[i] = unsigned_to_signed(x[i]);
        call_maxReduce<int>(tmp, n, comm_rank_of_max);
        for (int i = 0; i < n; i++)
            x[i] = signed_to_unsigned(tmp[i]);
        delete[] tmp;
    }
}
// int
template <>
void MPI_CLASS::call_maxReduce<int>(const int *x, int *y, int n,
                                    int *comm_rank_of_max) const {
    PROFILE_START("maxReduce1<int>", profile_level);
    if (comm_rank_of_max == nullptr) {
        MPI_Allreduce((void *)x, (void *)y, n, MPI_INT, MPI_MAX, communicator);
    } else {
        auto recv = new IntIntStruct[n];
        auto send = new IntIntStruct[n];
        for (int i = 0; i < n; ++i) {
            send[i].j = x[i];
            send[i].i = comm_rank;
        }
        MPI_Allreduce(send, recv, n, MPI_2INT, MPI_MAXLOC, communicator);
        for (int i = 0; i < n; ++i) {
            y[i] = recv[i].j;
            comm_rank_of_max[i] = recv[i].i;
        }
        delete[] recv;
        delete[] send;
    }
    PROFILE_STOP("maxReduce1<int>", profile_level);
}
template <>
void MPI_CLASS::call_maxReduce<int>(int *x, int n,
                                    int *comm_rank_of_max) const {
    PROFILE_START("maxReduce2<int>", profile_level);
    if (comm_rank_of_max == nullptr) {
        int *send = x;
        auto recv = new int[n];
        MPI_Allreduce(send, recv, n, MPI_INT, MPI_MAX, communicator);
        for (int i = 0; i < n; i++)
            x[i] = recv[i];
        delete[] recv;
    } else {
        auto recv = new IntIntStruct[n];
        auto send = new IntIntStruct[n];
        for (int i = 0; i < n; ++i) {
            send[i].j = x[i];
            send[i].i = comm_rank;
        }
        MPI_Allreduce(send, recv, n, MPI_2INT, MPI_MAXLOC, communicator);
        for (int i = 0; i < n; ++i) {
            x[i] = recv[i].j;
            comm_rank_of_max[i] = recv[i].i;
        }
        delete[] recv;
        delete[] send;
    }
    PROFILE_STOP("maxReduce2<int>", profile_level);
}
// long int
template <>
void MPI_CLASS::call_maxReduce<long int>(const long int *x, long int *y, int n,
                                         int *comm_rank_of_max) const {
    PROFILE_START("maxReduce1<lond int>", profile_level);
    if (comm_rank_of_max == nullptr) {
        MPI_Allreduce((void *)x, (void *)y, n, MPI_LONG, MPI_MAX, communicator);
    } else {
        auto recv = new LongIntStruct[n];
        auto send = new LongIntStruct[n];
        for (int i = 0; i < n; ++i) {
            send[i].j = x[i];
            send[i].i = comm_rank;
        }
        MPI_Allreduce(send, recv, n, MPI_LONG_INT, MPI_MAXLOC, communicator);
        for (int i = 0; i < n; ++i) {
            y[i] = recv[i].j;
            comm_rank_of_max[i] = recv[i].i;
        }
        delete[] recv;
        delete[] send;
    }
    PROFILE_STOP("maxReduce1<lond int>", profile_level);
}
template <>
void MPI_CLASS::call_maxReduce<long int>(long int *x, int n,
                                         int *comm_rank_of_max) const {
    PROFILE_START("maxReduce2<lond int>", profile_level);
    if (comm_rank_of_max == nullptr) {
        auto send = x;
        auto recv = new long int[n];
        MPI_Allreduce(send, recv, n, MPI_LONG, MPI_MAX, communicator);
        for (int i = 0; i < n; i++)
            x[i] = recv[i];
        delete[] recv;
    } else {
        auto recv = new LongIntStruct[n];
        auto send = new LongIntStruct[n];
        for (int i = 0; i < n; ++i) {
            send[i].j = x[i];
            send[i].i = comm_rank;
        }
        MPI_Allreduce(send, recv, n, MPI_LONG_INT, MPI_MAXLOC, communicator);
        for (int i = 0; i < n; ++i) {
            x[i] = recv[i].j;
            comm_rank_of_max[i] = recv[i].i;
        }
        delete[] recv;
        delete[] send;
    }
    PROFILE_STOP("maxReduce2<lond int>", profile_level);
}
// unsigned long int
template <>
void MPI_CLASS::call_maxReduce<unsigned long int>(const unsigned long int *send,
                                                  unsigned long int *recv,
                                                  int n,
                                                  int *comm_rank_of_max) const {
    if (comm_rank_of_max == nullptr) {
        PROFILE_START("maxReduce1<unsigned long>", profile_level);
        MPI_Allreduce((void *)send, (void *)recv, n, MPI_UNSIGNED_LONG, MPI_MAX,
                      communicator);
        PROFILE_STOP("maxReduce1<unsigned long>", profile_level);
    } else {
        auto tmp = new long int[n];
        for (int i = 0; i < n; i++)
            tmp[i] = unsigned_to_signed(send[i]);
        call_maxReduce<long int>(tmp, n, comm_rank_of_max);
        for (int i = 0; i < n; i++)
            recv[i] = signed_to_unsigned(tmp[i]);
        delete[] tmp;
    }
}
template <>
void MPI_CLASS::call_maxReduce<unsigned long int>(unsigned long int *x, int n,
                                                  int *comm_rank_of_max) const {
    if (comm_rank_of_max == nullptr) {
        PROFILE_START("maxReduce2<unsigned long>", profile_level);
        auto send = x;
        auto recv = new unsigned long int[n];
        MPI_Allreduce(send, recv, n, MPI_UNSIGNED_LONG, MPI_MAX, communicator);
        for (int i = 0; i < n; i++)
            x[i] = recv[i];
        delete[] recv;
        PROFILE_STOP("maxReduce2<unsigned long>", profile_level);
    } else {
        auto tmp = new long int[n];
        for (int i = 0; i < n; i++)
            tmp[i] = unsigned_to_signed(x[i]);
        call_maxReduce<long int>(tmp, n, comm_rank_of_max);
        for (int i = 0; i < n; i++)
            x[i] = signed_to_unsigned(tmp[i]);
        delete[] tmp;
    }
}
// unsigned long long int
template <>
void MPI_CLASS::call_maxReduce<unsigned long long int>(
    const unsigned long long int *send, unsigned long long int *recv, int n,
    int *comm_rank_of_max) const {
    PROFILE_START("maxReduce1<long int>", profile_level);
    if (comm_rank_of_max == nullptr) {
        auto x = new long long int[n];
        auto y = new long long int[n];
        for (int i = 0; i < n; i++)
            x[i] = unsigned_to_signed(send[i]);
        MPI_Allreduce((void *)x, (void *)y, n, MPI_LONG_LONG_INT, MPI_MAX,
                      communicator);
        for (int i = 0; i < n; i++)
            recv[i] = signed_to_unsigned(y[i]);
        delete[] x;
        delete[] y;
    } else {
        printf("maxReduce<long long int> will use double\n");
        auto tmp = new double[n];
        for (int i = 0; i < n; i++)
            tmp[i] = static_cast<double>(send[i]);
        call_maxReduce<double>(tmp, n, comm_rank_of_max);
        for (int i = 0; i < n; i++)
            recv[i] = static_cast<long long int>(tmp[i]);
        delete[] tmp;
    }
    PROFILE_STOP("maxReduce1<long int>", profile_level);
}
template <>
void MPI_CLASS::call_maxReduce<unsigned long long int>(
    unsigned long long int *x, int n, int *comm_rank_of_max) const {
    auto recv = new unsigned long long int[n];
    call_maxReduce<unsigned long long int>(x, recv, n, comm_rank_of_max);
    for (int i = 0; i < n; i++)
        x[i] = recv[i];
    delete[] recv;
}
// long long int
template <>
void MPI_CLASS::call_maxReduce<long long int>(const long long int *x,
                                              long long int *y, int n,
                                              int *comm_rank_of_max) const {
    PROFILE_START("maxReduce1<long int>", profile_level);
    if (comm_rank_of_max == nullptr) {
        MPI_Allreduce((void *)x, (void *)y, n, MPI_LONG_LONG_INT, MPI_MAX,
                      communicator);
    } else {
        printf("maxReduce<long long int> will use double\n");
        auto tmp = new double[n];
        for (int i = 0; i < n; i++)
            tmp[i] = static_cast<double>(x[i]);
        call_maxReduce<double>(tmp, n, comm_rank_of_max);
        for (int i = 0; i < n; i++)
            y[i] = static_cast<long long int>(tmp[i]);
        delete[] tmp;
    }
    PROFILE_STOP("maxReduce1<long int>", profile_level);
}
template <>
void MPI_CLASS::call_maxReduce<long long int>(long long int *x, int n,
                                              int *comm_rank_of_max) const {
    auto recv = new long long int[n];
    call_maxReduce<long long int>(x, recv, n, comm_rank_of_max);
    for (int i = 0; i < n; i++)
        x[i] = signed_to_unsigned(recv[i]);
    delete[] recv;
}
// float
template <>
void MPI_CLASS::call_maxReduce<float>(const float *x, float *y, int n,
                                      int *comm_rank_of_max) const {
    PROFILE_START("maxReduce1<float>", profile_level);
    if (comm_rank_of_max == nullptr) {
        MPI_Allreduce((void *)x, (void *)y, n, MPI_FLOAT, MPI_MAX,
                      communicator);
    } else {
        auto recv = new FloatIntStruct[n];
        auto send = new FloatIntStruct[n];
        for (int i = 0; i < n; ++i) {
            send[i].f = x[i];
            send[i].i = comm_rank;
        }
        MPI_Allreduce(send, recv, n, MPI_FLOAT_INT, MPI_MAXLOC, communicator);
        for (int i = 0; i < n; ++i) {
            y[i] = recv[i].f;
            comm_rank_of_max[i] = recv[i].i;
        }
        delete[] recv;
        delete[] send;
    }
    PROFILE_STOP("maxReduce1<float>", profile_level);
}
template <>
void MPI_CLASS::call_maxReduce<float>(float *x, int n,
                                      int *comm_rank_of_max) const {
    PROFILE_START("maxReduce2<float>", profile_level);
    if (comm_rank_of_max == nullptr) {
        auto send = x;
        auto recv = new float[n];
        MPI_Allreduce(send, recv, n, MPI_FLOAT, MPI_MAX, communicator);
        for (int i = 0; i < n; i++)
            x[i] = recv[i];
        delete[] recv;
    } else {
        auto recv = new FloatIntStruct[n];
        auto send = new FloatIntStruct[n];
        for (int i = 0; i < n; ++i) {
            send[i].f = x[i];
            send[i].i = comm_rank;
        }
        MPI_Allreduce(send, recv, n, MPI_FLOAT_INT, MPI_MAXLOC, communicator);
        for (int i = 0; i < n; ++i) {
            x[i] = recv[i].f;
            comm_rank_of_max[i] = recv[i].i;
        }
        delete[] recv;
        delete[] send;
    }
    PROFILE_STOP("maxReduce2<float>", profile_level);
}
// double
template <>
void MPI_CLASS::call_maxReduce<double>(const double *x, double *y, int n,
                                       int *comm_rank_of_max) const {
    PROFILE_START("maxReduce1<double>", profile_level);
    if (comm_rank_of_max == nullptr) {
        MPI_Allreduce((void *)x, (void *)y, n, MPI_DOUBLE, MPI_MAX,
                      communicator);
    } else {
        auto recv = new DoubleIntStruct[n];
        auto send = new DoubleIntStruct[n];
        for (int i = 0; i < n; ++i) {
            send[i].d = x[i];
            send[i].i = comm_rank;
        }
        MPI_Allreduce(send, recv, n, MPI_DOUBLE_INT, MPI_MAXLOC, communicator);
        for (int i = 0; i < n; ++i) {
            y[i] = recv[i].d;
            comm_rank_of_max[i] = recv[i].i;
        }
        delete[] recv;
        delete[] send;
    }
    PROFILE_STOP("maxReduce1<double>", profile_level);
}
template <>
void MPI_CLASS::call_maxReduce<double>(double *x, int n,
                                       int *comm_rank_of_max) const {
    PROFILE_START("maxReduce2<double>", profile_level);
    if (comm_rank_of_max == nullptr) {
        auto send = x;
        auto recv = new double[n];
        MPI_Allreduce(send, recv, n, MPI_DOUBLE, MPI_MAX, communicator);
        for (int i = 0; i < n; i++)
            x[i] = recv[i];
        delete[] recv;
    } else {
        auto recv = new DoubleIntStruct[n];
        auto send = new DoubleIntStruct[n];
        for (int i = 0; i < n; ++i) {
            send[i].d = x[i];
            send[i].i = comm_rank;
        }
        MPI_Allreduce(send, recv, n, MPI_DOUBLE_INT, MPI_MAXLOC, communicator);
        for (int i = 0; i < n; ++i) {
            x[i] = recv[i].d;
            comm_rank_of_max[i] = recv[i].i;
        }
        delete[] recv;
        delete[] send;
    }
    PROFILE_STOP("maxReduce2<double>", profile_level);
}
#endif

/************************************************************************
 *  bcast                                                                *
 *  Note: these specializations are only called when using MPI.          *
 ************************************************************************/
#ifdef USE_MPI
// char
template <>
void MPI_CLASS::call_bcast<unsigned char>(unsigned char *x, int n,
                                          int root) const {
    PROFILE_START("bcast<unsigned char>", profile_level);
    MPI_Bcast(x, n, MPI_UNSIGNED_CHAR, root, communicator);
    PROFILE_STOP("bcast<unsigned char>", profile_level);
}
template <> void MPI_CLASS::call_bcast<char>(char *x, int n, int root) const {
    PROFILE_START("bcast<char>", profile_level);
    MPI_Bcast(x, n, MPI_CHAR, root, communicator);
    PROFILE_STOP("bcast<char>", profile_level);
}
// int
template <>
void MPI_CLASS::call_bcast<unsigned int>(unsigned int *x, int n,
                                         int root) const {
    PROFILE_START("bcast<unsigned int>", profile_level);
    MPI_Bcast(x, n, MPI_UNSIGNED, root, communicator);
    PROFILE_STOP("bcast<unsigned int>", profile_level);
}
template <> void MPI_CLASS::call_bcast<int>(int *x, int n, int root) const {
    PROFILE_START("bcast<int>", profile_level);
    MPI_Bcast(x, n, MPI_INT, root, communicator);
    PROFILE_STOP("bcast<int>", profile_level);
}
// float
template <> void MPI_CLASS::call_bcast<float>(float *x, int n, int root) const {
    PROFILE_START("bcast<float>", profile_level);
    MPI_Bcast(x, n, MPI_FLOAT, root, communicator);
    PROFILE_STOP("bcast<float>", profile_level);
}
// double
template <>
void MPI_CLASS::call_bcast<double>(double *x, int n, int root) const {
    PROFILE_START("bcast<double>", profile_level);
    MPI_Bcast(x, n, MPI_DOUBLE, root, communicator);
    PROFILE_STOP("bcast<double>", profile_level);
}
#else
// We need a concrete instantiation of bcast<char>(x,n,root);
template <> void MPI_CLASS::call_bcast<char>(char *, int, int) const {}
#endif

/************************************************************************
 *  Perform a global barrier across all processors.                      *
 ************************************************************************/
void MPI_CLASS::barrier() const {
#ifdef USE_MPI
    MPI_Barrier(communicator);
#endif
}

/************************************************************************
 *  Send data array to another processor.                                *
 *  Note: these specializations are only called when using MPI.          *
 ************************************************************************/
#ifdef USE_MPI
// char
template <>
void MPI_CLASS::send<char>(const char *buf, int length, int recv_proc_number,
                           int tag) const {
    // Set the tag to 0 if it is < 0
    tag = (tag >= 0) ? tag : 0;
    MPI_INSIST(tag <= d_maxTag, "Maximum tag value exceeded");
    // Send the data
    PROFILE_START("send<char>", profile_level);
    MPI_Send((void *)buf, length, MPI_CHAR, recv_proc_number, tag,
             communicator);
    PROFILE_STOP("send<char>", profile_level);
}
// int
template <>
void MPI_CLASS::send<int>(const int *buf, int length, int recv_proc_number,
                          int tag) const {
    // Set the tag to 0 if it is < 0
    tag = (tag >= 0) ? tag : 0;
    MPI_INSIST(tag <= d_maxTag, "Maximum tag value exceeded");
    // Send the data
    PROFILE_START("send<int>", profile_level);
    MPI_Send((void *)buf, length, MPI_INT, recv_proc_number, tag, communicator);
    PROFILE_STOP("send<int>", profile_level);
}
// float
template <>
void MPI_CLASS::send<float>(const float *buf, int length, int recv_proc_number,
                            int tag) const {
    // Set the tag to 0 if it is < 0
    tag = (tag >= 0) ? tag : 0;
    MPI_INSIST(tag <= d_maxTag, "Maximum tag value exceeded");
    // Send the data
    PROFILE_START("send<float>", profile_level);
    MPI_Send((void *)buf, length, MPI_FLOAT, recv_proc_number, tag,
             communicator);
    PROFILE_STOP("send<float>", profile_level);
}
// double
template <>
void MPI_CLASS::send<double>(const double *buf, int length,
                             int recv_proc_number, int tag) const {
    // Set the tag to 0 if it is < 0
    tag = (tag >= 0) ? tag : 0;
    MPI_INSIST(tag <= d_maxTag, "Maximum tag value exceeded");
    // Send the data
    PROFILE_START("send<double>", profile_level);
    MPI_Send((void *)buf, length, MPI_DOUBLE, recv_proc_number, tag,
             communicator);
    PROFILE_STOP("send<double>", profile_level);
}
#else
// We need a concrete instantiation of send for use without MPI
template <>
void MPI_CLASS::send<char>(const char *buf, int length, int, int tag) const {
    MPI_INSIST(tag <= d_maxTag, "Maximum tag value exceeded");
    MPI_INSIST(tag >= 0, "tag must be >= 0");
    PROFILE_START("send<char>", profile_level);
    auto id = getRequest(communicator, tag);
    auto it = global_isendrecv_list.find(id);
    MPI_INSIST(it == global_isendrecv_list.end(),
               "send must be paired with a previous call to irecv in serial");
    MPI_ASSERT(it->second.status == 2);
    memcpy((char *)it->second.data, buf, length);
    global_isendrecv_list.erase(it);
    PROFILE_START("send<char>", profile_level);
}
#endif

/************************************************************************
 *  Non-blocking send data array to another processor.                   *
 *  Note: these specializations are only called when using MPI.          *
 ************************************************************************/
#ifdef USE_MPI
// char
template <>
MPI_Request MPI_CLASS::Isend<char>(const char *buf, int length, int recv_proc,
                                   int tag) const {
    MPI_INSIST(tag <= d_maxTag, "Maximum tag value exceeded");
    MPI_INSIST(tag >= 0, "tag must be >= 0");
    MPI_Request request;
    PROFILE_START("Isend<char>", profile_level);
    MPI_Isend((void *)buf, length, MPI_CHAR, recv_proc, tag, communicator,
              &request);
    PROFILE_STOP("Isend<char>", profile_level);
    return request;
}
// int
template <>
MPI_Request MPI_CLASS::Isend<int>(const int *buf, int length, int recv_proc,
                                  int tag) const {
    MPI_INSIST(tag <= d_maxTag, "Maximum tag value exceeded");
    MPI_INSIST(tag >= 0, "tag must be >= 0");
    MPI_Request request;
    PROFILE_START("Isend<int>", profile_level);
    MPI_Isend((void *)buf, length, MPI_INT, recv_proc, tag, communicator,
              &request);
    PROFILE_STOP("Isend<int>", profile_level);
    return request;
}
// float
template <>
MPI_Request MPI_CLASS::Isend<float>(const float *buf, int length, int recv_proc,
                                    int tag) const {
    MPI_INSIST(tag <= d_maxTag, "Maximum tag value exceeded");
    MPI_INSIST(tag >= 0, "tag must be >= 0");
    MPI_Request request;
    PROFILE_START("Isend<float>", profile_level);
    MPI_Isend((void *)buf, length, MPI_FLOAT, recv_proc, tag, communicator,
              &request);
    PROFILE_STOP("Isend<float>", profile_level);
    return request;
}
// double
template <>
MPI_Request MPI_CLASS::Isend<double>(const double *buf, int length,
                                     int recv_proc, int tag) const {
    MPI_INSIST(tag <= d_maxTag, "Maximum tag value exceeded");
    MPI_INSIST(tag >= 0, "tag must be >= 0");
    MPI_Request request;
    PROFILE_START("Isend<double>", profile_level);
    MPI_Isend((void *)buf, length, MPI_DOUBLE, recv_proc, tag, communicator,
              &request);
    PROFILE_STOP("Isend<double>", profile_level);
    return request;
}
#else
// We need a concrete instantiation of send for use without mpi
template <>
MPI_Request MPI_CLASS::Isend<char>(const char *buf, int length, int,
                                   int tag) const {
    MPI_INSIST(tag <= d_maxTag, "Maximum tag value exceeded");
    MPI_INSIST(tag >= 0, "tag must be >= 0");
    PROFILE_START("Isend<char>", profile_level);
    auto id = getRequest(communicator, tag);
    auto it = global_isendrecv_list.find(id);
    if (it == global_isendrecv_list.end()) {
        // We are calling isend first
        Isendrecv_struct data;
        data.data = buf;
        data.status = 1;
        global_isendrecv_list.insert(
            std::pair<MPI_Request, Isendrecv_struct>(id, data));
    } else {
        // We called irecv first
        MPI_ASSERT(it->second.status == 2);
        memcpy((char *)it->second.data, buf, length);
        global_isendrecv_list.erase(it);
    }
    PROFILE_STOP("Isend<char>", profile_level);
    return id;
}
#endif

/************************************************************************
 *  Send byte array to another processor.                                *
 ************************************************************************/
void MPI_CLASS::sendBytes(const void *buf, int number_bytes,
                          int recv_proc_number, int tag) const {
    MPI_INSIST(tag <= d_maxTag, "Maximum tag value exceeded");
    MPI_INSIST(tag >= 0, "tag must be >= 0");
    send<char>((const char *)buf, number_bytes, recv_proc_number, tag);
}

/************************************************************************
 *  Non-blocking send byte array to another processor.                   *
 ************************************************************************/
MPI_Request MPI_CLASS::IsendBytes(const void *buf, int number_bytes,
                                  const int recv_proc, const int tag) const {
    MPI_INSIST(tag <= d_maxTag, "Maximum tag value exceeded");
    MPI_INSIST(tag >= 0, "tag must be >= 0");
    return Isend<char>((const char *)buf, number_bytes, recv_proc, tag);
}

/************************************************************************
 *  Recieve data array to another processor.                             *
 *  Note: these specializations are only called when using MPI.          *
 ************************************************************************/
#ifdef USE_MPI
// char
template <>
void MPI_CLASS::recv<char>(char *buf, int &length, int send_proc_number,
                           const bool get_length, int tag) const {
    // Set the tag to 0 if it is < 0
    tag = (tag >= 0) ? tag : 0;
    MPI_INSIST(tag <= d_maxTag, "Maximum tag value exceeded");
    PROFILE_START("recv<char>", profile_level);
    // Get the recieve length if necessary
    if (get_length) {
        int bytes = this->probe(send_proc_number, tag);
        int recv_length = bytes / sizeof(char);
        MPI_INSIST(length >= recv_length,
                   "Recived length is larger than allocated array");
        length = recv_length;
    }
    // Send the data
    MPI_Status status;
    MPI_Recv((void *)buf, length, MPI_CHAR, send_proc_number, tag, communicator,
             &status);
    PROFILE_STOP("recv<char>", profile_level);
}
// int
template <>
void MPI_CLASS::recv<int>(int *buf, int &length, int send_proc_number,
                          const bool get_length, int tag) const {
    // Set the tag to 0 if it is < 0
    tag = (tag >= 0) ? tag : 0;
    MPI_INSIST(tag <= d_maxTag, "Maximum tag value exceeded");
    PROFILE_START("recv<int>", profile_level);
    // Get the recieve length if necessary
    if (get_length) {
        int bytes = this->probe(send_proc_number, tag);
        int recv_length = bytes / sizeof(int);
        MPI_INSIST(length >= recv_length,
                   "Recived length is larger than allocated array");
        length = recv_length;
    }
    // Send the data
    MPI_Status status;
    MPI_Recv((void *)buf, length, MPI_INT, send_proc_number, tag, communicator,
             &status);
    PROFILE_STOP("recv<int>", profile_level);
}
// float
template <>
void MPI_CLASS::recv<float>(float *buf, int &length, int send_proc_number,
                            const bool get_length, int tag) const {
    // Set the tag to 0 if it is < 0
    tag = (tag >= 0) ? tag : 0;
    MPI_INSIST(tag <= d_maxTag, "Maximum tag value exceeded");
    PROFILE_START("recv<float>", profile_level);
    // Get the recieve length if necessary
    if (get_length) {
        int bytes = this->probe(send_proc_number, tag);
        int recv_length = bytes / sizeof(float);
        MPI_INSIST(length >= recv_length,
                   "Recived length is larger than allocated array");
        length = recv_length;
    }
    // Send the data
    MPI_Status status;
    MPI_Recv((void *)buf, length, MPI_FLOAT, send_proc_number, tag,
             communicator, &status);
    PROFILE_STOP("recv<float>", profile_level);
}
// double
template <>
void MPI_CLASS::recv<double>(double *buf, int &length, int send_proc_number,
                             const bool get_length, int tag) const {
    // Set the tag to 0 if it is < 0
    tag = (tag >= 0) ? tag : 0;
    MPI_INSIST(tag <= d_maxTag, "Maximum tag value exceeded");
    PROFILE_START("recv<double>", profile_level);
    // Get the recieve length if necessary
    if (get_length) {
        int bytes = this->probe(send_proc_number, tag);
        int recv_length = bytes / sizeof(double);
        MPI_INSIST(length >= recv_length,
                   "Recived length is larger than allocated array");
        length = recv_length;
    }
    // Send the data
    MPI_Status status;
    MPI_Recv((void *)buf, length, MPI_DOUBLE, send_proc_number, tag,
             communicator, &status);
    PROFILE_STOP("recv<double>", profile_level);
}
#else
// We need a concrete instantiation of recv for use without mpi
template <>
void MPI_CLASS::recv<char>(char *buf, int &length, int, const bool,
                           int tag) const {
    MPI_INSIST(tag <= d_maxTag, "Maximum tag value exceeded");
    MPI_INSIST(tag >= 0, "tag must be >= 0");
    PROFILE_START("recv<char>", profile_level);
    auto id = getRequest(communicator, tag);
    auto it = global_isendrecv_list.find(id);
    MPI_INSIST(it != global_isendrecv_list.end(),
               "recv must be paired with a previous call to isend in serial");
    MPI_ASSERT(it->second.status == 1);
    memcpy(buf, it->second.data, length);
    global_isendrecv_list.erase(it);
    PROFILE_STOP("recv<char>", profile_level);
}
#endif

/************************************************************************
 *  Non-blocking recieve data array to another processor.                *
 *  Note: these specializations are only called when using MPI.          *
 ************************************************************************/
#ifdef USE_MPI
// char
template <>
MPI_Request MPI_CLASS::Irecv<char>(char *buf, int length, int send_proc,
                                   int tag) const {
    MPI_INSIST(tag <= d_maxTag, "Maximum tag value exceeded");
    MPI_INSIST(tag >= 0, "tag must be >= 0");
    MPI_Request request;
    PROFILE_START("Irecv<char>", profile_level);
    MPI_Irecv((void *)buf, length, MPI_CHAR, send_proc, tag, communicator,
              &request);
    PROFILE_STOP("Irecv<char>", profile_level);
    return request;
}
// int
template <>
MPI_Request MPI_CLASS::Irecv<int>(int *buf, int length, int send_proc,
                                  int tag) const {
    MPI_INSIST(tag <= d_maxTag, "Maximum tag value exceeded");
    MPI_INSIST(tag >= 0, "tag must be >= 0");
    MPI_Request request;
    PROFILE_START("Irecv<int>", profile_level);
    MPI_Irecv((void *)buf, length, MPI_INT, send_proc, tag, communicator,
              &request);
    PROFILE_STOP("Irecv<int>", profile_level);
    return request;
}
// float
template <>
MPI_Request MPI_CLASS::Irecv<float>(float *buf, int length, int send_proc,
                                    int tag) const {
    MPI_INSIST(tag <= d_maxTag, "Maximum tag value exceeded");
    MPI_INSIST(tag >= 0, "tag must be >= 0");
    MPI_Request request;
    PROFILE_START("Irecv<float>", profile_level);
    MPI_Irecv((void *)buf, length, MPI_FLOAT, send_proc, tag, communicator,
              &request);
    PROFILE_STOP("Irecv<float>", profile_level);
    return request;
}
// double
template <>
MPI_Request MPI_CLASS::Irecv<double>(double *buf, int length, int send_proc,
                                     int tag) const {
    MPI_INSIST(tag <= d_maxTag, "Maximum tag value exceeded");
    MPI_INSIST(tag >= 0, "tag must be >= 0");
    MPI_Request request;
    PROFILE_START("Irecv<double>", profile_level);
    MPI_Irecv((void *)buf, length, MPI_DOUBLE, send_proc, tag, communicator,
              &request);
    PROFILE_STOP("Irecv<double>", profile_level);
    return request;
}
#else
// We need a concrete instantiation of irecv for use without mpi
template <>
MPI_Request MPI_CLASS::Irecv<char>(char *buf, int length, int, int tag) const {
    MPI_INSIST(tag <= d_maxTag, "Maximum tag value exceeded");
    MPI_INSIST(tag >= 0, "tag must be >= 0");
    PROFILE_START("Irecv<char>", profile_level);
    auto id = getRequest(communicator, tag);
    auto it = global_isendrecv_list.find(id);
    if (it == global_isendrecv_list.end()) {
        // We are calling Irecv first
        Isendrecv_struct data;
        data.data = buf;
        data.status = 2;
        global_isendrecv_list.insert(
            std::pair<MPI_Request, Isendrecv_struct>(id, data));
    } else {
        // We called Isend first
        MPI_ASSERT(it->second.status == 1);
        memcpy(buf, it->second.data, length);
        global_isendrecv_list.erase(it);
    }
    PROFILE_STOP("Irecv<char>", profile_level);
    return id;
}
#endif

/************************************************************************
 *  Recieve byte array to another processor.                             *
 ************************************************************************/
void MPI_CLASS::recvBytes(void *buf, int &number_bytes, int send_proc,
                          int tag) const {
    recv<char>((char *)buf, number_bytes, send_proc, false, tag);
}

/************************************************************************
 *  Recieve byte array to another processor.                             *
 ************************************************************************/
MPI_Request MPI_CLASS::IrecvBytes(void *buf, int number_bytes, int send_proc,
                                  int tag) const {
    MPI_INSIST(tag <= d_maxTag, "Maximum tag value exceeded");
    MPI_INSIST(tag >= 0, "tag must be >= 0");
    return Irecv<char>((char *)buf, number_bytes, send_proc, tag);
}

/************************************************************************
 *  sendrecv                                                             *
 ************************************************************************/
#if defined(USE_MPI)
template <>
void MPI_CLASS::sendrecv<char>(const char *sendbuf, int sendcount, int dest,
                               int sendtag, char *recvbuf, int recvcount,
                               int source, int recvtag) const {
    PROFILE_START("sendrecv<char>", profile_level);
    MPI_Sendrecv(sendbuf, sendcount, MPI_CHAR, dest, sendtag, recvbuf,
                 recvcount, MPI_CHAR, source, recvtag, communicator,
                 MPI_STATUS_IGNORE);
    PROFILE_STOP("sendrecv<char>", profile_level);
}
template <>
void MPI_CLASS::sendrecv<int>(const int *sendbuf, int sendcount, int dest,
                              int sendtag, int *recvbuf, int recvcount,
                              int source, int recvtag) const {
    PROFILE_START("sendrecv<int>", profile_level);
    MPI_Sendrecv(sendbuf, sendcount, MPI_INT, dest, sendtag, recvbuf, recvcount,
                 MPI_INT, source, recvtag, communicator, MPI_STATUS_IGNORE);
    PROFILE_STOP("sendrecv<int>", profile_level);
}
template <>
void MPI_CLASS::sendrecv<float>(const float *sendbuf, int sendcount, int dest,
                                int sendtag, float *recvbuf, int recvcount,
                                int source, int recvtag) const {
    PROFILE_START("sendrecv<float>", profile_level);
    MPI_Sendrecv(sendbuf, sendcount, MPI_FLOAT, dest, sendtag, recvbuf,
                 recvcount, MPI_FLOAT, source, recvtag, communicator,
                 MPI_STATUS_IGNORE);
    PROFILE_STOP("sendrecv<float>", profile_level);
}
template <>
void MPI_CLASS::sendrecv<double>(const double *sendbuf, int sendcount, int dest,
                                 int sendtag, double *recvbuf, int recvcount,
                                 int source, int recvtag) const {
    PROFILE_START("sendrecv<double>", profile_level);
    MPI_Sendrecv(sendbuf, sendcount, MPI_DOUBLE, dest, sendtag, recvbuf,
                 recvcount, MPI_DOUBLE, source, recvtag, communicator,
                 MPI_STATUS_IGNORE);
    PROFILE_STOP("sendrecv<double>", profile_level);
}
#endif

/************************************************************************
 *  allGather                                                            *
 *  Note: these specializations are only called when using MPI.          *
 ************************************************************************/
#ifdef USE_MPI
// unsigned char
template <>
void MPI_CLASS::call_allGather<unsigned char>(const unsigned char &x_in,
                                              unsigned char *x_out) const {
    PROFILE_START("allGather<unsigned char>", profile_level);
    MPI_Allgather((void *)&x_in, 1, MPI_UNSIGNED_CHAR, (void *)x_out, 1,
                  MPI_UNSIGNED_CHAR, communicator);
    PROFILE_STOP("allGather<unsigned char>", profile_level);
}
template <>
void MPI_CLASS::call_allGather<unsigned char>(const unsigned char *x_in,
                                              int size_in, unsigned char *x_out,
                                              int *size_out,
                                              int *disp_out) const {
    PROFILE_START("allGatherv<unsigned char>", profile_level);
    MPI_Allgatherv((void *)x_in, size_in, MPI_CHAR, (void *)x_out, size_out,
                   disp_out, MPI_CHAR, communicator);
    PROFILE_STOP("allGatherv<unsigned char>", profile_level);
}
// char
template <>
void MPI_CLASS::call_allGather<char>(const char &x_in, char *x_out) const {
    PROFILE_START("allGather<char>", profile_level);
    MPI_Allgather((void *)&x_in, 1, MPI_CHAR, (void *)x_out, 1, MPI_CHAR,
                  communicator);
    PROFILE_STOP("allGather<char>", profile_level);
}
template <>
void MPI_CLASS::call_allGather<char>(const char *x_in, int size_in, char *x_out,
                                     int *size_out, int *disp_out) const {
    PROFILE_START("allGatherv<char>", profile_level);
    MPI_Allgatherv((void *)x_in, size_in, MPI_CHAR, (void *)x_out, size_out,
                   disp_out, MPI_CHAR, communicator);
    PROFILE_STOP("allGatherv<char>", profile_level);
}
// unsigned int
template <>
void MPI_CLASS::call_allGather<unsigned int>(const unsigned int &x_in,
                                             unsigned int *x_out) const {
    PROFILE_START("allGather<unsigned int>", profile_level);
    MPI_Allgather((void *)&x_in, 1, MPI_UNSIGNED, (void *)x_out, 1,
                  MPI_UNSIGNED, communicator);
    PROFILE_STOP("allGather<unsigned int>", profile_level);
}
template <>
void MPI_CLASS::call_allGather<unsigned int>(const unsigned int *x_in,
                                             int size_in, unsigned int *x_out,
                                             int *size_out,
                                             int *disp_out) const {
    PROFILE_START("allGatherv<unsigned int>", profile_level);
    MPI_Allgatherv((void *)x_in, size_in, MPI_UNSIGNED, (void *)x_out, size_out,
                   disp_out, MPI_UNSIGNED, communicator);
    PROFILE_STOP("allGatherv<unsigned int>", profile_level);
}
// int
template <>
void MPI_CLASS::call_allGather<int>(const int &x_in, int *x_out) const {
    PROFILE_START("allGather<int>", profile_level);
    MPI_Allgather((void *)&x_in, 1, MPI_INT, (void *)x_out, 1, MPI_INT,
                  communicator);
    PROFILE_STOP("allGather<int>", profile_level);
}
template <>
void MPI_CLASS::call_allGather<int>(const int *x_in, int size_in, int *x_out,
                                    int *size_out, int *disp_out) const {
    PROFILE_START("allGatherv<int>", profile_level);
    MPI_Allgatherv((void *)x_in, size_in, MPI_INT, (void *)x_out, size_out,
                   disp_out, MPI_INT, communicator);
    PROFILE_STOP("allGatherv<int>", profile_level);
}
// unsigned long int
template <>
void MPI_CLASS::call_allGather<unsigned long int>(
    const unsigned long int &x_in, unsigned long int *x_out) const {
    PROFILE_START("allGather<unsigned long>", profile_level);
    MPI_Allgather((void *)&x_in, 1, MPI_UNSIGNED_LONG, (void *)x_out, 1,
                  MPI_UNSIGNED_LONG, communicator);
    PROFILE_STOP("allGather<unsigned long>", profile_level);
}
template <>
void MPI_CLASS::call_allGather<unsigned long int>(const unsigned long int *x_in,
                                                  int size_in,
                                                  unsigned long int *x_out,
                                                  int *size_out,
                                                  int *disp_out) const {
    PROFILE_START("allGatherv<unsigned long>", profile_level);
    MPI_Allgatherv((void *)x_in, size_in, MPI_UNSIGNED_LONG, (void *)x_out,
                   size_out, disp_out, MPI_UNSIGNED_LONG, communicator);
    PROFILE_STOP("allGatherv<unsigned long>", profile_level);
}
// long int
template <>
void MPI_CLASS::call_allGather<long int>(const long int &x_in,
                                         long int *x_out) const {
    PROFILE_START("allGather<long int>", profile_level);
    MPI_Allgather((void *)&x_in, 1, MPI_LONG, (void *)x_out, 1, MPI_LONG,
                  communicator);
    PROFILE_STOP("allGather<long int>", profile_level);
}
template <>
void MPI_CLASS::call_allGather<long int>(const long int *x_in, int size_in,
                                         long int *x_out, int *size_out,
                                         int *disp_out) const {
    PROFILE_START("allGatherv<long int>", profile_level);
    MPI_Allgatherv((void *)x_in, size_in, MPI_LONG, (void *)x_out, size_out,
                   disp_out, MPI_LONG, communicator);
    PROFILE_STOP("allGatherv<long int>", profile_level);
}
// float
template <>
void MPI_CLASS::call_allGather<float>(const float &x_in, float *x_out) const {
    PROFILE_START("allGather<float>", profile_level);
    MPI_Allgather((void *)&x_in, 1, MPI_FLOAT, (void *)x_out, 1, MPI_FLOAT,
                  communicator);
    PROFILE_STOP("allGather<float>", profile_level);
}
template <>
void MPI_CLASS::call_allGather<float>(const float *x_in, int size_in,
                                      float *x_out, int *size_out,
                                      int *disp_out) const {
    PROFILE_START("allGatherv<float>", profile_level);
    MPI_Allgatherv((void *)x_in, size_in, MPI_FLOAT, (void *)x_out, size_out,
                   disp_out, MPI_FLOAT, communicator);
    PROFILE_STOP("allGatherv<float>", profile_level);
}
// double
template <>
void MPI_CLASS::call_allGather<double>(const double &x_in,
                                       double *x_out) const {
    PROFILE_START("allGather<double>", profile_level);
    MPI_Allgather((void *)&x_in, 1, MPI_DOUBLE, (void *)x_out, 1, MPI_DOUBLE,
                  communicator);
    PROFILE_STOP("allGather<double>", profile_level);
}
template <>
void MPI_CLASS::call_allGather<double>(const double *x_in, int size_in,
                                       double *x_out, int *size_out,
                                       int *disp_out) const {
    PROFILE_START("allGatherv<double>", profile_level);
    MPI_Allgatherv((void *)x_in, size_in, MPI_DOUBLE, (void *)x_out, size_out,
                   disp_out, MPI_DOUBLE, communicator);
    PROFILE_STOP("allGatherv<double>", profile_level);
}
#else
// We need a concrete instantiation of call_allGather<char>(x_in,size_in,x_out,size_out)
template <>
void MPI_CLASS::call_allGather<char>(const char *, int, char *, int *,
                                     int *) const {
    MPI_ERROR("Internal error in communicator (allGather) ");
}
#endif

/************************************************************************
 *  allToAll                                                             *
 *  Note: these specializations are only called when using MPI.          *
 ************************************************************************/
#ifdef USE_MPI
template <>
void MPI_CLASS::allToAll<unsigned char>(int n, const unsigned char *send,
                                        unsigned char *recv) const {
    PROFILE_START("allToAll<unsigned char>", profile_level);
    MPI_Alltoall((void *)send, n, MPI_UNSIGNED_CHAR, (void *)recv, n,
                 MPI_UNSIGNED_CHAR, communicator);
    PROFILE_STOP("allToAll<unsigned char>", profile_level);
}
template <>
void MPI_CLASS::allToAll<char>(int n, const char *send, char *recv) const {
    PROFILE_START("allToAll<char>", profile_level);
    MPI_Alltoall((void *)send, n, MPI_CHAR, (void *)recv, n, MPI_CHAR,
                 communicator);
    PROFILE_STOP("allToAll<char>", profile_level);
}
template <>
void MPI_CLASS::allToAll<unsigned int>(int n, const unsigned int *send,
                                       unsigned int *recv) const {
    PROFILE_START("allToAll<unsigned int>", profile_level);
    MPI_Alltoall((void *)send, n, MPI_UNSIGNED, (void *)recv, n, MPI_UNSIGNED,
                 communicator);
    PROFILE_STOP("allToAll<unsigned int>", profile_level);
}
template <>
void MPI_CLASS::allToAll<int>(int n, const int *send, int *recv) const {
    PROFILE_START("allToAll<int>", profile_level);
    MPI_Alltoall((void *)send, n, MPI_INT, (void *)recv, n, MPI_INT,
                 communicator);
    PROFILE_STOP("allToAll<int>", profile_level);
}
template <>
void MPI_CLASS::allToAll<unsigned long int>(int n,
                                            const unsigned long int *send,
                                            unsigned long int *recv) const {
    PROFILE_START("allToAll<unsigned long>", profile_level);
    MPI_Alltoall((void *)send, n, MPI_UNSIGNED_LONG, (void *)recv, n,
                 MPI_UNSIGNED_LONG, communicator);
    PROFILE_STOP("allToAll<unsigned long>", profile_level);
}
template <>
void MPI_CLASS::allToAll<long int>(int n, const long int *send,
                                   long int *recv) const {
    PROFILE_START("allToAll<long int>", profile_level);
    MPI_Alltoall((void *)send, n, MPI_LONG, (void *)recv, n, MPI_LONG,
                 communicator);
    PROFILE_STOP("allToAll<long int>", profile_level);
}
template <>
void MPI_CLASS::allToAll<float>(int n, const float *send, float *recv) const {
    PROFILE_START("allToAll<float>", profile_level);
    MPI_Alltoall((void *)send, n, MPI_FLOAT, (void *)recv, n, MPI_FLOAT,
                 communicator);
    PROFILE_STOP("allToAll<float>", profile_level);
}
template <>
void MPI_CLASS::allToAll<double>(int n, const double *send,
                                 double *recv) const {
    PROFILE_START("allToAll<double>", profile_level);
    MPI_Alltoall((void *)send, n, MPI_DOUBLE, (void *)recv, n, MPI_DOUBLE,
                 communicator);
    PROFILE_STOP("allToAll<double>", profile_level);
}
#endif

/************************************************************************
 *  call_allToAll                                                        *
 *  Note: these specializations are only called when using MPI.          *
 ************************************************************************/
#ifdef USE_MPI
// unsigned char
template <>
void MPI_CLASS::call_allToAll<unsigned char>(
    const unsigned char *send_data, const int send_cnt[], const int send_disp[],
    unsigned char *recv_data, const int *recv_cnt, const int *recv_disp) const {
    PROFILE_START("allToAllv<unsigned char>", profile_level);
    MPI_Alltoallv((void *)send_data, (int *)send_cnt, (int *)send_disp,
                  MPI_UNSIGNED_CHAR, (void *)recv_data, (int *)recv_cnt,
                  (int *)recv_disp, MPI_UNSIGNED_CHAR, communicator);
    PROFILE_STOP("allToAllv<unsigned char>", profile_level);
}
// char
template <>
void MPI_CLASS::call_allToAll<char>(const char *send_data, const int send_cnt[],
                                    const int send_disp[], char *recv_data,
                                    const int *recv_cnt,
                                    const int *recv_disp) const {
    PROFILE_START("allToAllv<char>", profile_level);
    MPI_Alltoallv((void *)send_data, (int *)send_cnt, (int *)send_disp,
                  MPI_CHAR, (void *)recv_data, (int *)recv_cnt,
                  (int *)recv_disp, MPI_CHAR, communicator);
    PROFILE_STOP("allToAllv<char>", profile_level);
}
// unsigned int
template <>
void MPI_CLASS::call_allToAll<unsigned int>(
    const unsigned int *send_data, const int send_cnt[], const int send_disp[],
    unsigned int *recv_data, const int *recv_cnt, const int *recv_disp) const {
    PROFILE_START("allToAllv<unsigned int>", profile_level);
    MPI_Alltoallv((void *)send_data, (int *)send_cnt, (int *)send_disp,
                  MPI_UNSIGNED, (void *)recv_data, (int *)recv_cnt,
                  (int *)recv_disp, MPI_UNSIGNED, communicator);
    PROFILE_STOP("allToAllv<unsigned int>", profile_level);
}
// int
template <>
void MPI_CLASS::call_allToAll<int>(const int *send_data, const int send_cnt[],
                                   const int send_disp[], int *recv_data,
                                   const int *recv_cnt,
                                   const int *recv_disp) const {
    PROFILE_START("allToAllv<int>", profile_level);
    MPI_Alltoallv((void *)send_data, (int *)send_cnt, (int *)send_disp, MPI_INT,
                  (void *)recv_data, (int *)recv_cnt, (int *)recv_disp, MPI_INT,
                  communicator);
    PROFILE_STOP("allToAllv<int>", profile_level);
}
// unsigned long int
template <>
void MPI_CLASS::call_allToAll<unsigned long int>(
    const unsigned long int *send_data, const int send_cnt[],
    const int send_disp[], unsigned long int *recv_data, const int *recv_cnt,
    const int *recv_disp) const {
    PROFILE_START("allToAllv<unsigned long>", profile_level);
    MPI_Alltoallv((void *)send_data, (int *)send_cnt, (int *)send_disp,
                  MPI_UNSIGNED_LONG, (void *)recv_data, (int *)recv_cnt,
                  (int *)recv_disp, MPI_UNSIGNED_LONG, communicator);
    PROFILE_STOP("allToAllv<unsigned long>", profile_level);
}
// long int
template <>
void MPI_CLASS::call_allToAll<long int>(
    const long int *send_data, const int send_cnt[], const int send_disp[],
    long int *recv_data, const int *recv_cnt, const int *recv_disp) const {
    PROFILE_START("allToAllv<long int>", profile_level);
    MPI_Alltoallv((void *)send_data, (int *)send_cnt, (int *)send_disp,
                  MPI_LONG, (void *)recv_data, (int *)recv_cnt,
                  (int *)recv_disp, MPI_LONG, communicator);
    PROFILE_STOP("allToAllv<long int>", profile_level);
}
// float
template <>
void MPI_CLASS::call_allToAll<float>(const float *send_data,
                                     const int send_cnt[],
                                     const int send_disp[], float *recv_data,
                                     const int *recv_cnt,
                                     const int *recv_disp) const {
    PROFILE_START("allToAllv<float>", profile_level);
    MPI_Alltoallv((void *)send_data, (int *)send_cnt, (int *)send_disp,
                  MPI_FLOAT, (void *)recv_data, (int *)recv_cnt,
                  (int *)recv_disp, MPI_FLOAT, communicator);
    PROFILE_STOP("allToAllv<float>", profile_level);
}
// double
template <>
void MPI_CLASS::call_allToAll<double>(const double *send_data,
                                      const int send_cnt[],
                                      const int send_disp[], double *recv_data,
                                      const int *recv_cnt,
                                      const int *recv_disp) const {
    PROFILE_START("allToAllv<double>", profile_level);
    MPI_Alltoallv((void *)send_data, (int *)send_cnt, (int *)send_disp,
                  MPI_DOUBLE, (void *)recv_data, (int *)recv_cnt,
                  (int *)recv_disp, MPI_DOUBLE, communicator);
    PROFILE_STOP("allToAllv<double>", profile_level);
}
#else
// Default instatiation of unsigned char
template <>
void MPI_CLASS::call_allToAll<char>(const char *, const int[], const int[],
                                    char *, const int *, const int *) const {
    MPI_ERROR("Should not reach this point");
}
#endif

/************************************************************************
 *  call_sumScan                                                         *
 *  Note: these specializations are only called when using MPI.          *
 ************************************************************************/
#ifdef USE_MPI
// unsigned char
template <>
void MPI_CLASS::call_sumScan<unsigned char>(const unsigned char *send,
                                            unsigned char *recv, int n) const {
    PROFILE_START("sumScan<unsigned char>", profile_level);
    MPI_Scan((void *)send, (void *)recv, n, MPI_UNSIGNED_CHAR, MPI_SUM,
             communicator);
    PROFILE_STOP("sumScan<unsigned char>", profile_level);
}
// char
template <>
void MPI_CLASS::call_sumScan<char>(const char *send, char *recv, int n) const {
    PROFILE_START("sumScan<char>", profile_level);
    MPI_Scan((void *)send, (void *)recv, n, MPI_SIGNED_CHAR, MPI_SUM,
             communicator);
    PROFILE_STOP("sumScan<char>", profile_level);
}
// unsigned int
template <>
void MPI_CLASS::call_sumScan<unsigned int>(const unsigned int *send,
                                           unsigned int *recv, int n) const {
    PROFILE_START("sumScan<unsigned int>", profile_level);
    MPI_Scan((void *)send, (void *)recv, n, MPI_UNSIGNED, MPI_SUM,
             communicator);
    PROFILE_STOP("sumScan<unsigned int>", profile_level);
}
// int
template <>
void MPI_CLASS::call_sumScan<int>(const int *send, int *recv, int n) const {
    PROFILE_START("sumScan<int>", profile_level);
    MPI_Scan((void *)send, (void *)recv, n, MPI_INT, MPI_SUM, communicator);
    PROFILE_STOP("sumScan<int>", profile_level);
}
// long int
template <>
void MPI_CLASS::call_sumScan<long int>(const long int *send, long int *recv,
                                       int n) const {
    PROFILE_START("sumScan<long int>", profile_level);
    MPI_Scan((void *)send, (void *)recv, n, MPI_LONG, MPI_SUM, communicator);
    PROFILE_STOP("sumScan<long int>", profile_level);
}
// unsigned long int
template <>
void MPI_CLASS::call_sumScan<unsigned long>(const unsigned long *send,
                                            unsigned long *recv, int n) const {
    PROFILE_START("sumScan<unsigned long>", profile_level);
    MPI_Scan((void *)send, (void *)recv, n, MPI_UNSIGNED_LONG, MPI_SUM,
             communicator);
    PROFILE_STOP("sumScan<unsigned long>", profile_level);
}
// size_t
#ifdef USE_WINDOWS
template <>
void MPI_CLASS::call_sumScan<size_t>(const size_t *send, size_t *recv,
                                     int n) const {
    MPI_ASSERT(MPI_SIZE_T != 0);
    PROFILE_START("sumScan<size_t>", profile_level);
    MPI_Scan((void *)send, (void *)recv, n, MPI_SIZE_T, MPI_SUM, communicator);
    PROFILE_STOP("sumScan<size_t>", profile_level);
}
#endif
// float
template <>
void MPI_CLASS::call_sumScan<float>(const float *send, float *recv,
                                    int n) const {
    PROFILE_START("sumScan<float>", profile_level);
    MPI_Scan((void *)send, (void *)recv, n, MPI_FLOAT, MPI_SUM, communicator);
    PROFILE_STOP("sumScan<float>", profile_level);
}
// double
template <>
void MPI_CLASS::call_sumScan<double>(const double *send, double *recv,
                                     int n) const {
    PROFILE_START("sumScan<double>", profile_level);
    MPI_Scan((void *)send, (void *)recv, n, MPI_DOUBLE, MPI_SUM, communicator);
    PROFILE_STOP("sumScan<double>", profile_level);
}
// std::complex<double>
template <>
void MPI_CLASS::call_sumScan<std::complex<double>>(
    const std::complex<double> *x, std::complex<double> *y, int n) const {
    auto send = new double[2 * n];
    auto recv = new double[2 * n];
    for (int i = 0; i < n; i++) {
        send[2 * i + 0] = real(x[i]);
        send[2 * i + 1] = imag(x[i]);
    }
    MPI_Scan((void *)send, (void *)recv, 2 * n, MPI_DOUBLE, MPI_SUM,
             communicator);
    for (int i = 0; i < n; i++)
        y[i] = std::complex<double>(recv[2 * i + 0], recv[2 * i + 1]);
    delete[] send;
    delete[] recv;
}
#endif

/************************************************************************
 *  call_minScan                                                         *
 *  Note: these specializations are only called when using MPI.          *
 ************************************************************************/
#ifdef USE_MPI
// unsigned char
template <>
void MPI_CLASS::call_minScan<unsigned char>(const unsigned char *send,
                                            unsigned char *recv, int n) const {
    PROFILE_START("minScan<unsigned char>", profile_level);
    MPI_Scan((void *)send, (void *)recv, n, MPI_UNSIGNED_CHAR, MPI_MIN,
             communicator);
    PROFILE_STOP("minScan<unsigned char>", profile_level);
}
// char
template <>
void MPI_CLASS::call_minScan<char>(const char *send, char *recv, int n) const {
    PROFILE_START("minScan<char>", profile_level);
    MPI_Scan((void *)send, (void *)recv, n, MPI_SIGNED_CHAR, MPI_MIN,
             communicator);
    PROFILE_STOP("minScan<char>", profile_level);
}
// unsigned int
template <>
void MPI_CLASS::call_minScan<unsigned int>(const unsigned int *send,
                                           unsigned int *recv, int n) const {
    PROFILE_START("minScan<unsigned int>", profile_level);
    MPI_Scan((void *)send, (void *)recv, n, MPI_UNSIGNED, MPI_MIN,
             communicator);
    PROFILE_STOP("minScan<unsigned int>", profile_level);
}
// int
template <>
void MPI_CLASS::call_minScan<int>(const int *send, int *recv, int n) const {
    PROFILE_START("minScan<int>", profile_level);
    MPI_Scan((void *)send, (void *)recv, n, MPI_INT, MPI_MIN, communicator);
    PROFILE_STOP("minScan<int>", profile_level);
}
// unsigned long int
template <>
void MPI_CLASS::call_minScan<unsigned long int>(const unsigned long int *send,
                                                unsigned long int *recv,
                                                int n) const {
    PROFILE_START("minScan<unsigned long>", profile_level);
    MPI_Scan((void *)send, (void *)recv, n, MPI_UNSIGNED_LONG, MPI_MIN,
             communicator);
    PROFILE_STOP("minScan<unsigned long>", profile_level);
}
// long int
template <>
void MPI_CLASS::call_minScan<long int>(const long int *send, long int *recv,
                                       int n) const {
    PROFILE_START("minScan<long int>", profile_level);
    MPI_Scan((void *)send, (void *)recv, n, MPI_LONG, MPI_MIN, communicator);
    PROFILE_STOP("minScan<long int>", profile_level);
}
// size_t
#ifdef USE_WINDOWS
template <>
void MPI_CLASS::call_minScan<size_t>(const size_t *send, size_t *recv,
                                     int n) const {
    MPI_ASSERT(MPI_SIZE_T != 0);
    PROFILE_START("minScan<size_t>", profile_level);
    MPI_Scan((void *)send, (void *)recv, n, MPI_SIZE_T, MPI_MIN, communicator);
    PROFILE_STOP("minScan<size_t>", profile_level);
}
#endif
// float
template <>
void MPI_CLASS::call_minScan<float>(const float *send, float *recv,
                                    int n) const {
    PROFILE_START("minScan<float>", profile_level);
    MPI_Scan((void *)send, (void *)recv, n, MPI_FLOAT, MPI_MIN, communicator);
    PROFILE_STOP("minScan<float>", profile_level);
}
// double
template <>
void MPI_CLASS::call_minScan<double>(const double *send, double *recv,
                                     int n) const {
    PROFILE_START("minScan<double>", profile_level);
    MPI_Scan((void *)send, (void *)recv, n, MPI_DOUBLE, MPI_MIN, communicator);
    PROFILE_STOP("minScan<double>", profile_level);
}
#endif

/************************************************************************
 *  call_maxScan                                                         *
 *  Note: these specializations are only called when using MPI.          *
 ************************************************************************/
#ifdef USE_MPI
// unsigned char
template <>
void MPI_CLASS::call_maxScan<unsigned char>(const unsigned char *send,
                                            unsigned char *recv, int n) const {
    PROFILE_START("maxScan<unsigned char>", profile_level);
    MPI_Scan((void *)send, (void *)recv, n, MPI_UNSIGNED_CHAR, MPI_MAX,
             communicator);
    PROFILE_STOP("maxScan<unsigned char>", profile_level);
}
// char
template <>
void MPI_CLASS::call_maxScan<char>(const char *send, char *recv, int n) const {
    PROFILE_START("maxScan<char>", profile_level);
    MPI_Scan((void *)send, (void *)recv, n, MPI_SIGNED_CHAR, MPI_MAX,
             communicator);
    PROFILE_STOP("maxScan<char>", profile_level);
}
// unsigned int
template <>
void MPI_CLASS::call_maxScan<unsigned int>(const unsigned int *send,
                                           unsigned int *recv, int n) const {
    PROFILE_START("maxScan<unsigned int>", profile_level);
    MPI_Scan((void *)send, (void *)recv, n, MPI_UNSIGNED, MPI_MAX,
             communicator);
    PROFILE_STOP("maxScan<unsigned int>", profile_level);
}
// int
template <>
void MPI_CLASS::call_maxScan<int>(const int *send, int *recv, int n) const {
    PROFILE_START("maxScan<int>", profile_level);
    MPI_Scan((void *)send, (void *)recv, n, MPI_INT, MPI_MAX, communicator);
    PROFILE_STOP("maxScan<int>", profile_level);
}
// long int
template <>
void MPI_CLASS::call_maxScan<long int>(const long int *send, long int *recv,
                                       int n) const {
    PROFILE_START("maxScan<long int>", profile_level);
    MPI_Scan((void *)send, (void *)recv, n, MPI_LONG, MPI_MAX, communicator);
    PROFILE_STOP("maxScan<long int>", profile_level);
}
// unsigned long int
template <>
void MPI_CLASS::call_maxScan<unsigned long int>(const unsigned long int *send,
                                                unsigned long int *recv,
                                                int n) const {
    PROFILE_START("maxScan<unsigned long>", profile_level);
    MPI_Scan((void *)send, (void *)recv, n, MPI_UNSIGNED_LONG, MPI_MAX,
             communicator);
    PROFILE_STOP("maxScan<unsigned long>", profile_level);
}
// size_t
#ifdef USE_WINDOWS
template <>
void MPI_CLASS::call_maxScan<size_t>(const size_t *send, size_t *recv,
                                     int n) const {
    MPI_ASSERT(MPI_SIZE_T != 0);
    PROFILE_START("maxScan<size_t>", profile_level);
    MPI_Scan((void *)send, (void *)recv, n, MPI_SIZE_T, MPI_MAX, communicator);
    PROFILE_STOP("maxScan<size_t>", profile_level);
}
#endif
// float
template <>
void MPI_CLASS::call_maxScan<float>(const float *send, float *recv,
                                    int n) const {
    PROFILE_START("maxScan<float>", profile_level);
    MPI_Scan((void *)send, (void *)recv, n, MPI_INT, MPI_MAX, communicator);
    PROFILE_STOP("maxScan<float>", profile_level);
}
// double
template <>
void MPI_CLASS::call_maxScan<double>(const double *send, double *recv,
                                     int n) const {
    PROFILE_START("maxScan<double>", profile_level);
    MPI_Scan((void *)send, (void *)recv, n, MPI_DOUBLE, MPI_MAX, communicator);
    PROFILE_STOP("maxScan<double>", profile_level);
}
#endif

/************************************************************************
 *  Communicate ranks for communication                                  *
 ************************************************************************/
std::vector<int> MPI_CLASS::commRanks(const std::vector<int> &ranks) const {
#ifdef USE_MPI
    // Get a byte array with the ranks to communicate
    auto data1 = new char[comm_size];
    auto data2 = new char[comm_size];
    memset(data1, 0, comm_size);
    memset(data2, 0, comm_size);
    for (auto &rank : ranks)
        data1[rank] = 1;
    MPI_Alltoall(data1, 1, MPI_CHAR, data2, 1, MPI_CHAR, communicator);
    int N = 0;
    for (int i = 0; i < comm_size; i++)
        N += data2[i];
    std::vector<int> ranks_out;
    ranks_out.reserve(N);
    for (int i = 0; i < comm_size; i++) {
        if (data2[i])
            ranks_out.push_back(i);
    }
    delete[] data1;
    delete[] data2;
    return ranks_out;
#else
    return ranks;
#endif
}

/************************************************************************
 *  Wait functions                                                       *
 ************************************************************************/
#ifdef USE_MPI
void MPI_CLASS::wait(MPI_Request request) {
    PROFILE_START("wait", profile_level);
    MPI_Status status;
    MPI_Wait(&request, &status);
    /*int flag = 0;
    int err  = MPI_Test( &request, &flag, &status );
    MPI_ASSERT( err == MPI_SUCCESS ); // Check that the first call is valid
    while ( !flag ) {
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
        // Check if the request has finished
        MPI_Test( &request, &flag, &status );
    }*/
    PROFILE_STOP("wait", profile_level);
}
int MPI_CLASS::waitAny(int count, MPI_Request *request) {
    if (count == 0)
        return -1;
    PROFILE_START("waitAny", profile_level);
    int index = -1;
    auto status = new MPI_Status[count];
    MPI_Waitany(count, request, &index, status);
    /*int flag    = 0;
    int err     = MPI_Testany( count, request, &index, &flag, status );
    MPI_ASSERT( err == MPI_SUCCESS ); // Check that the first call is valid
    while ( !flag ) {
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
        // Check if the request has finished
        MPI_Testany( count, request, &index, &flag, status );
    }
    MPI_ASSERT( index >= 0 ); // Check that the index is valid*/
    delete[] status;
    PROFILE_STOP("waitAny", profile_level);
    return index;
}
void MPI_CLASS::waitAll(int count, MPI_Request *request) {
    if (count == 0)
        return;
    PROFILE_START("waitAll", profile_level);
    auto status = new MPI_Status[count];
    MPI_Waitall(count, request, status);
    /*int flag    = 0;
    int err     = MPI_Testall( count, request, &flag, status );
    MPI_ASSERT( err == MPI_SUCCESS ); // Check that the first call is valid
    while ( !flag ) {
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
        // Check if the request has finished
        MPI_Testall( count, request, &flag, status );
    }*/
    PROFILE_STOP("waitAll", profile_level);
    delete[] status;
}
std::vector<int> MPI_CLASS::waitSome(int count, MPI_Request *request) {
    if (count == 0)
        return std::vector<int>();
    PROFILE_START("waitSome", profile_level);
    std::vector<int> indicies(count, -1);
    auto *status = new MPI_Status[count];
    int outcount = 0;
    MPI_Waitsome(count, request, &outcount, indicies.data(), status);
    /*int err      = MPI_Testsome( count, request, &outcount, &indicies[0], status );
    MPI_ASSERT( err == MPI_SUCCESS );        // Check that the first call is valid
    MPI_ASSERT( outcount != MPI_UNDEFINED ); // Check that the first call is valid
    while ( outcount == 0 ) {
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
        // Check if the request has finished
        MPI_Testsome( count, request, &outcount, &indicies[0], status );
    }*/
    indicies.resize(outcount);
    delete[] status;
    PROFILE_STOP("waitSome", profile_level);
    return indicies;
}
#else
void MPI_CLASS::wait(MPI_Request request) {
    PROFILE_START("wait", profile_level);
    while (1) {
        // Check if the request is in our list
        if (global_isendrecv_list.find(request) == global_isendrecv_list.end())
            break;
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
    }
    PROFILE_STOP("wait", profile_level);
}
int MPI_CLASS::waitAny(int count, MPI_Request *request) {
    if (count == 0)
        return -1;
    PROFILE_START("waitAny", profile_level);
    int index = 0;
    while (1) {
        // Check if the request is in our list
        bool found_any = false;
        for (int i = 0; i < count; i++) {
            if (global_isendrecv_list.find(request[i]) ==
                global_isendrecv_list.end()) {
                found_any = true;
                index = i;
            }
        }
        if (found_any)
            break;
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
    }
    PROFILE_STOP("waitAny", profile_level);
    return index;
}
void MPI_CLASS::waitAll(int count, MPI_Request *request) {
    if (count == 0)
        return;
    PROFILE_START("waitAll", profile_level);
    while (1) {
        // Check if the request is in our list
        bool found_all = true;
        for (int i = 0; i < count; i++) {
            if (global_isendrecv_list.find(request[i]) !=
                global_isendrecv_list.end())
                found_all = false;
        }
        if (found_all)
            break;
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
    }
    PROFILE_STOP("waitAll", profile_level);
}
std::vector<int> MPI_CLASS::waitSome(int count, MPI_Request *request) {
    if (count == 0)
        return std::vector<int>();
    PROFILE_START("waitSome", profile_level);
    std::vector<int> indicies;
    while (1) {
        // Check if the request is in our list
        for (int i = 0; i < count; i++) {
            if (global_isendrecv_list.find(request[i]) ==
                global_isendrecv_list.end())
                indicies.push_back(i);
        }
        if (!indicies.empty())
            break;
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
    }
    PROFILE_STOP("waitSome", profile_level);
    return indicies;
}
#endif

/************************************************************************
 *  Probe functions                                                      *
 ************************************************************************/
#ifdef USE_MPI
int MPI_CLASS::Iprobe(int source, int tag) const {
    MPI_INSIST(tag <= d_maxTag, "Maximum tag value exceeded");
    MPI_INSIST(tag >= 0, "tag must be >= 0");
    MPI_Status status;
    int flag = 0;
    MPI_Iprobe(source, tag, communicator, &flag, &status);
    if (flag == 0)
        return -1;
    int count;
    MPI_Get_count(&status, MPI_BYTE, &count);
    MPI_ASSERT(count >= 0);
    return count;
}
int MPI_CLASS::probe(int source, int tag) const {
    MPI_INSIST(tag <= d_maxTag, "Maximum tag value exceeded");
    MPI_INSIST(tag >= 0, "tag must be >= 0");
    MPI_Status status;
    MPI_Probe(source, tag, communicator, &status);
    int count;
    MPI_Get_count(&status, MPI_BYTE, &count);
    MPI_ASSERT(count >= 0);
    return count;
}
#else
int MPI_CLASS::Iprobe(int, int) const {
    MPI_ERROR("Not implimented for serial codes (Iprobe)");
    return 0;
}
int MPI_CLASS::probe(int, int) const {
    MPI_ERROR("Not implimented for serial codes (probe)");
    return 0;
}
#endif

/************************************************************************
 *  Timer functions                                                      *
 ************************************************************************/
#ifdef USE_MPI
double MPI_CLASS::time() { return MPI_Wtime(); }
double MPI_CLASS::tick() { return MPI_Wtick(); }
#else
double MPI_CLASS::time() {
    auto t = std::chrono::system_clock::now();
    auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
        t.time_since_epoch());
    return 1e-9 * ns.count();
}
double MPI_CLASS::tick() {
    auto period = std::chrono::system_clock::period();
    return static_cast<double>(period.num) / static_cast<double>(period.den);
}
#endif

/************************************************************************
 *  Serialize a block of code across MPI processes                       *
 ************************************************************************/
void MPI_CLASS::serializeStart() {
#ifdef USE_MPI
    using namespace std::chrono_literals;
    if (comm_rank == 0) {
        // Start rank 0 immediately
    } else {
        // Wait for a message from the previous rank
        MPI_Request request;
        MPI_Status status;
        int flag = false, buf = 0;
        MPI_Irecv(&buf, 1, MPI_INT, comm_rank - 1, 5627, MPI_COMM_WORLD,
                  &request);
        while (!flag) {
            MPI_Test(&request, &flag, &status);
            std::this_thread::sleep_for(50ms);
        }
    }
#endif
}
void MPI_CLASS::serializeStop() {
#ifdef USE_MPI
    using namespace std::chrono_literals;
    if (comm_rank < comm_size - 1) {
        // Send flag to next rank
        MPI_Send(&comm_rank, 1, MPI_INT, comm_rank + 1, 5627, MPI_COMM_WORLD);
        // Wait for final finished flag
        int flag = false, buf = 0;
        MPI_Request request;
        MPI_Status status;
        MPI_Irecv(&buf, 1, MPI_INT, comm_size - 1, 5627, MPI_COMM_WORLD,
                  &request);
        while (!flag) {
            MPI_Test(&request, &flag, &status);
            std::this_thread::sleep_for(50ms);
        }
    } else {
        // Send final flag to all ranks
        for (int i = 0; i < comm_size - 1; i++)
            MPI_Send(&comm_rank, 1, MPI_INT, i, 5627, MPI_COMM_WORLD);
    }
#endif
}

/****************************************************************************
 * Function to start/stop MPI                                                *
 ****************************************************************************/
#ifdef USE_MPI
static bool called_MPI_Init = false;
#endif
bool MPI_CLASS::MPI_Active() {
#ifdef USE_MPI
    int MPI_initialized, MPI_finialized;
    MPI_Initialized(&MPI_initialized);
    MPI_Finalized(&MPI_finialized);
    return MPI_initialized != 0 && MPI_finialized == 0;
#else
    return false;
#endif
}
void MPI_CLASS::start_MPI(int argc, char *argv[], int profile_level) {
    changeProfileLevel(profile_level);
    NULL_USE(argc);
    NULL_USE(argv);
#ifdef USE_MPI
    if (MPI_Active()) {
        called_MPI_Init = false;
    } else {
        int provided;
        int result =
            MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
        if (result != MPI_SUCCESS)
            MPI_ERROR("Unable to initialize MPI");
        if (provided < MPI_THREAD_MULTIPLE)
            std::cerr
                << "Warning: Failed to start MPI with MPI_THREAD_MULTIPLE\n";
        called_MPI_Init = true;
    }
#endif
}
void MPI_CLASS::stop_MPI() {
#ifdef USE_MPI
    int finalized;
    MPI_Finalized(&finalized);
    if (called_MPI_Init && !finalized) {
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        called_MPI_Init = true;
    }
#endif
}

/****************************************************************************
 * Function to perform load balancing                                        *
 ****************************************************************************/
MPI MPI::loadBalance(double local, std::vector<double> work) {
    MPI_ASSERT((int)work.size() == getSize());
    auto perf = allGather(local);
    std::vector<int> I(work.size());
    for (size_t i = 0; i < work.size(); i++)
        I[i] = i;
    auto J = I;
    quicksort(perf, I);
    quicksort(work, J);
    std::vector<int> key(work.size());
    for (size_t i = 0; i < work.size(); i++)
        key[J[i]] = I[i];
    return split(0, key[getRank()]);
}

/****************************************************************************
 * Function Persistent Communication                                         *
 ****************************************************************************/
template <>
std::shared_ptr<MPI_Request> MPI::Isend_init<double>(const double *buf, int N,
                                                     int proc, int tag) const {
    std::shared_ptr<MPI_Request> obj(new MPI_Request, [](MPI_Request *req) {
        MPI_Request_free(req);
        delete req;
    });
    MPI_Send_init(buf, N, MPI_DOUBLE, proc, tag, communicator, obj.get());
    return obj;
}
template <>
std::shared_ptr<MPI_Request> MPI::Irecv_init<double>(double *buf, int N,
                                                     int proc, int tag) const {
    std::shared_ptr<MPI_Request> obj(new MPI_Request, [](MPI_Request *req) {
        MPI_Request_free(req);
        delete req;
    });
    MPI_Recv_init(buf, N, MPI_DOUBLE, proc, tag, communicator, obj.get());
    return obj;
}
void MPI::Start(MPI_Request &request) { MPI_Start(&request); }

} // namespace Utilities

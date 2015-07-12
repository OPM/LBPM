#ifndef COMMON_H_INC
#define COMMON_H_INC

#include "common/Array.h"
#include "common/Communication.h"


/*!
 * @brief  Compute the blob
 * @details  Compute the blob (F>vf|S>vs) starting from (i,j,k) - oil blob
 * @return  Returns the number of cubes in the blob
 * @param[out] blobs        blobs
 * @param[out] nblobs       Number of blobs
 * @param[out] ncubes       Number of cubes
 * @param[out] indicator    indicator
 * @param[in]  F            F
 * @param[in]  S            S
 * @param[in]  vf           vf
 * @param[in]  vs           vs
 * @param[in]  startx       startx
 * @param[in]  starty       starty
 * @param[in/out] temp      temp
 */
int ComputeBlob( IntArray &blobs, int &nblobs, int &ncubes, IntArray &indicator,
	   const DoubleArray &F, const DoubleArray &S, double vf, double vs, int startx, int starty,
	   int startz, IntArray &temp, bool periodic=true );


/*!
 * @brief  Compute the blob
 * @details  Compute the blob (F>vf|S>vs) starting from (i,j,k) - oil blob
 * @return  Returns the number of cubes in the blob
 * @param[in] Phase         Phase
 * @param[in] SignDist      SignDist
 * @param[in] vF            vF
 * @param[in] vS            vS
 * @param[in] S             S
 * @param[out] LocalBlobID  The ids of the blobs
 * @return  Returns the number of blobs
 */
int ComputeLocalBlobIDs( const DoubleArray& Phase, const DoubleArray& SignDist, 
    double vF, double vS, IntArray& LocalBlobID, bool periodic=true );

/*!
 *  @brief Compute blob of an arbitrary phase
 *  @details Compute the connected components for Phase(i,j,k)=VALUE
 *  @return the number of connected components of the phase
 *  @param[in] PhaseID
 *  @param[in] VALUE
 *  @param[out] ComponentLabel
 *  @param[in] periodic
 */

int ComputeLocalPhaseComponent(IntArray &PhaseID, int VALUE, IntArray &ComponentLabel,
		bool periodic )

/*!
 * @brief  Compute the blob
 * @details  Compute the blob (F>vf|S>vs) starting from (i,j,k) - oil blob
 * @return  Returns the number of cubes in the blob
 * @param[in] nx            Number of elements in the x-direction
 * @param[in] ny            Number of elements in the y-direction
 * @param[in] nz            Number of elements in the z-direction
 * @param[in] Phase         Phase
 * @param[in] SignDist      SignDist
 * @param[in] vF            vF
 * @param[in] vS            vS
 * @param[in] S             S
 * @param[out] LocalBlobID  The ids of the blobs
 * @return  Returns the number of blobs
 */
int ComputeGlobalBlobIDs( int nx, int ny, int nz, RankInfoStruct rank_info, 
    const DoubleArray& Phase, const DoubleArray& SignDist, double vF, double vS, 
    IntArray& GlobalBlobID );


/*!
 * @brief  Reorder the blobs
 * @details  Reorder the blobs based on the number of cells they contain
 *    largest first.
 * @param[in] nx            Number of elements in the x-direction
 * @param[in] ny            Number of elements in the y-direction
 * @param[in] nz            Number of elements in the z-direction
 * @param[in/out] ID        The ids of the blobs
 */
void ReorderBlobIDs( IntArray& ID );



#endif

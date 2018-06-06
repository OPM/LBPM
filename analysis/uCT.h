#ifndef uCT_H_INC
#define uCT_H_INC

#include "common/Array.h"
#include "common/Domain.h"
#include "common/Communication.h"



/*!
 * @brief  Interpolate between meshes
 * @details  This routine interpolates from a coarse to a fine mesh
 * @param[in] Coarse    Coarse mesh solution
 * @param[out] Fine     Fine mesh solution
 */
void InterpolateMesh( const Array<float> &Coarse, Array<float> &Fine );


// Smooth the data using the distance
void smooth( const Array<float>& VOL, const Array<float>& Dist, float sigma, Array<float>& MultiScaleSmooth, fillHalo<float>& fillFloat );


// Segment the data
void segment( const Array<float>& data, Array<char>& ID, float tol );


// Remove disconnected phases
void removeDisconnected( Array<char>& ID, const Domain& Dm );


// Solve a level (without any coarse level information)
void solve( const Array<float>& VOL, Array<float>& Mean, Array<char>& ID,
    Array<float>& Dist, Array<float>& MultiScaleSmooth, Array<float>& NonLocalMean, 
    fillHalo<float>& fillFloat, const Domain& Dm, int nprocx,
    float threshold, float lamda, float sigsq, int depth);


// Refine a solution from a coarse grid to a fine grid
void refine( const Array<float>& Dist_coarse, 
    const Array<float>& VOL, Array<float>& Mean, Array<char>& ID,
    Array<float>& Dist, Array<float>& MultiScaleSmooth, Array<float>& NonLocalMean, 
    fillHalo<float>& fillFloat, const Domain& Dm, int nprocx, int level );


// Remove regions that are likely noise by shrinking the volumes by dx,
// removing all values that are more than dx+delta from the surface, and then
// growing by dx+delta and intersecting with the original data
void filter_final( Array<char>& ID, Array<float>& Dist,
    fillHalo<float>& fillFloat, const Domain& Dm,
    Array<float>& Mean, Array<float>& Dist1, Array<float>& Dist2 );


// Filter the original data
void filter_src( const Domain& Dm, Array<float>& src );


#endif

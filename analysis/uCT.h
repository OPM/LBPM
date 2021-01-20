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
    fillHalo<float>& fillFloat, const Domain& Dm, int nprocx, int level,
    float threshold, float lamda, float sigsq, int depth);


// Remove regions that are likely noise by shrinking the volumes by dx,
// removing all values that are more than dx+delta from the surface, and then
// growing by dx+delta and intersecting with the original data
void filter_final( Array<char>& ID, Array<float>& Dist,
    fillHalo<float>& fillFloat, const Domain& Dm,
    Array<float>& Mean, Array<float>& Dist1, Array<float>& Dist2 );


// Filter the original data
void filter_src( const Domain& Dm, Array<float>& src );


#endif

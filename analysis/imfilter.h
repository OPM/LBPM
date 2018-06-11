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
// These functions mimic the behavior of imfilter in MATLAB
#ifndef included_imfilter
#define included_imfilter

#include "common/Utilities.h"
#include "common/Array.h"
#include <vector>


namespace imfilter {


//! enum to store the BC type
enum class BC { fixed=0, symmetric=1, replicate=2, circular=3 };


/*!
 * @brief  N-D filtering of multidimensional images
 * @details  imfilter filters the multidimensional array A with the
 *    multidimensional filter H.  The result B has the same size and class as A.
 * @param[in]  A            The input array (Nx,Ny,Nz)
 * @param[in]  H            The filter (2*Nhx+1,2*Nhy+1,...)
 * @param[in]  boundary     The boundary conditions to apply (ndim):
 *                          fixed     - Input array values outside the bounds of the array are
 *                                      implicitly assumed to have the value X
 *                          symmetric - Input array values outside the bounds of the array are
 *                                      computed by mirror-reflecting the array across the array border
 *                          replicate - Input array values outside the bounds of the array are
 *                                      assumed to equal the nearest array border value
 *                          circular  - Input array values outside the bounds of the array are
 *                                      computed by implicitly assuming the input array is periodic.
 * @param[in]  X            The value to use for boundary conditions (only used if boundary==fixed)
 */
template<class TYPE>
Array<TYPE> imfilter( const Array<TYPE>& A, const Array<TYPE>& H, const std::vector<imfilter::BC>& boundary, const TYPE X=0 );


/*!
 * @brief  N-D filtering of multidimensional images
 * @details  imfilter filters the multidimensional array A with the
 *    multidimensional filter H.  The result B has the same size and class as A.
 * @param[in]  A            The input array (Nx,Ny,Nz)
 * @param[in]  Nh           The size of the filter
 * @param[in]  H            The filter function to use ( y = H(data) )
 *                          Note that the data passed to this function will be of size 2*Nh+1
 * @param[in]  boundary     The boundary conditions to apply (ndim):
 *                          fixed     - Input array values outside the bounds of the array are
 *                                      implicitly assumed to have the value X
 *                          symmetric - Input array values outside the bounds of the array are
 *                                      computed by mirror-reflecting the array across the array border
 *                          replicate - Input array values outside the bounds of the array are
 *                                      assumed to equal the nearest array border value
 *                          circular  - Input array values outside the bounds of the array are
 *                                      computed by implicitly assuming the input array is periodic.
 * @param[in]  X            The value to use for boundary conditions (only used if boundary==fixed)
 */
template<class TYPE>
Array<TYPE> imfilter( const Array<TYPE>& A, const std::vector<int>& Nh,
    std::function<TYPE(const Array<TYPE>&)> H,
    const std::vector<imfilter::BC>& boundary, const TYPE X=0 );


/*!
 * @brief  N-D filtering of multidimensional images
 * @details  imfilter filters the multidimensional array A with the
 *    multidimensional filter H.  The result B has the same size and class as A.
 *    This version works with separable filters and is more efficient than a single filter.
 * @param[in]  A            The input array (Nx,Ny,Nz)
 * @param[in]  H            The filter [2*Nhx+1,2*Nhy+1,...]
 * @param[in]  boundary     The boundary conditions to apply (ndim):
 *                          fixed     - Input array values outside the bounds of the array are
 *                                      implicitly assumed to have the value X
 *                          symmetric - Input array values outside the bounds of the array are
 *                                      computed by mirror-reflecting the array across the array border
 *                          replicate - Input array values outside the bounds of the array are
 *                                      assumed to equal the nearest array border value
 *                          circular  - Input array values outside the bounds of the array are
 *                                      computed by implicitly assuming the input array is periodic.
 * @param[in]  X            The value to use for boundary conditions (only used if boundary==fixed)
 */
template<class TYPE>
Array<TYPE> imfilter_separable( const Array<TYPE>& A, const std::vector<Array<TYPE>>& H,
    const std::vector<imfilter::BC>& boundary, const TYPE X=0 );


/*!
 * @brief  N-D filtering of multidimensional images
 * @details  imfilter filters the multidimensional array A with the
 *    multidimensional filter H.  The result B has the same size and class as A.
 *    This version works with separable filters and is more efficient than a single filter.
 * @param[in]  A            The input array (Nx,Ny,Nz)
 * @param[in]  H            The filter [2*Nhx+1,2*Nhy+1,...]
 * @param[in]  boundary     The boundary conditions to apply (ndim):
 *                          fixed     - Input array values outside the bounds of the array are
 *                                      implicitly assumed to have the value X
 *                          symmetric - Input array values outside the bounds of the array are
 *                                      computed by mirror-reflecting the array across the array border
 *                          replicate - Input array values outside the bounds of the array are
 *                                      assumed to equal the nearest array border value
 *                          circular  - Input array values outside the bounds of the array are
 *                                      computed by implicitly assuming the input array is periodic.
 * @param[in]  X            The value to use for boundary conditions (only used if boundary==fixed)
 */
template<class TYPE>
Array<TYPE> imfilter_separable( const Array<TYPE>& A, const std::vector<int>& Nh,
    std::vector<std::function<TYPE(const Array<TYPE>&)>> H,
    const std::vector<imfilter::BC>& boundary, const TYPE X=0 );


/*!
 * @brief  N-D filtering of multidimensional images
 * @details  imfilter filters the multidimensional array A with the
 *    multidimensional filter H.  The result B has the same size and class as A.
 *    This version works with separable filters and is more efficient than a single filter.
 * @param[in]  A            The input array (Nx,Ny,Nz)
 * @param[in]  H            The filter [2*Nhx+1,2*Nhy+1,...]
 * @param[in]  boundary     The boundary conditions to apply (ndim):
 *                          fixed     - Input array values outside the bounds of the array are
 *                                      implicitly assumed to have the value X
 *                          symmetric - Input array values outside the bounds of the array are
 *                                      computed by mirror-reflecting the array across the array border
 *                          replicate - Input array values outside the bounds of the array are
 *                                      assumed to equal the nearest array border value
 *                          circular  - Input array values outside the bounds of the array are
 *                                      computed by implicitly assuming the input array is periodic.
 * @param[in]  X            The value to use for boundary conditions (only used if boundary==fixed)
 */
template<class TYPE>
Array<TYPE> imfilter_separable( const Array<TYPE>& A, const std::vector<int>& Nh,
    std::vector<std::function<TYPE(int, const TYPE*)>> H,
    const std::vector<imfilter::BC>& boundary, const TYPE X=0 );


/**
 * @brief Create a filter to use with imfilter
 * @details  This function creates one of several predefined filters
 *   to use with imfilter.  The filter will always sum to 1.
 *   Note: this function allocates memory with the new command, the user must call delete.
 *
 * @param[in]  N            The stencil size in each direction
 * @param[in]  type         The type of filter to create
 *                          average  - Simple averaging filter
 *                          gaussian - Gaussian filter with given standard deviation.
 *                                     Optional argument is a double array of size ndim
 *                                     giving the standard deviation in each direction.
 *                                     A default value of 0.5 is used if not provided.
 * \param[in] args          An optional argument that some of the filters use
 */
template<class TYPE>
Array<TYPE> create_filter( const std::vector<int>& N, const std::string &type, const void *args = NULL );


}


#include "analysis/imfilter.hpp"

#endif


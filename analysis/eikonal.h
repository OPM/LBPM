#ifndef Eikonal_H_INC
#define Eikonal_H_INC

#include "common/Domain.h"


/*!
 * @brief  Calculate the distance solving the Eikonal equation
 * @details  This routine converts the data in the Distance array to a signed distance
 * by solving the equation df/dt = sign*(1-|grad f|), where Distance provides
 * the values of f on the mesh associated with domain Dm
 * It has been tested with segmented data initialized to values [-1,1]
 * and will converge toward the signed distance to the surface bounding the associated phases
 *
 * Reference:
 * Min C (2010) On reinitializing level set functions, Journal of Computational Physics	229
 *
 * @param[in/out] Distance  Distance function
 * @param[in] ID            Segmentation id
 * @param[in] DM            Domain information
 * @param[in] timesteps     Maximum number of timesteps to process
 * @return  Returns the global variation
 */
template<class TYPE>
TYPE Eikonal( Array<TYPE> &Distance, const Array<char> &ID, const Domain &Dm, int timesteps);


/*!
 * @brief  Calculate the distance using a simple method
 * @details  This routine calculates the distance using a very simple method working off the segmentation id.
 *
 * @param[in/out] Distance  Distance function
 * @param[in] ID            Segmentation id
 * @param[in] DM            Domain information
 * @return  Returns the global variation
 */
void CalcDist3D( Array<float> &Distance, const Array<char> &ID, const Domain &Dm );


/*!
 * @brief  Calculate the distance using a multi-level method
 * @details  This routine calculates the distance using a multi-grid method
 *
 * @return  Returns the number of cubes in the blob
 * @param[in/out] Distance  Distance function
 * @param[in] ID            Segmentation id
 * @param[in] DM            Domain information
 * @return  Returns the global variation
 */
void CalcDistMultiLevel( Array<float> &Distance, const Array<char> &ID, const Domain &Dm );


#endif

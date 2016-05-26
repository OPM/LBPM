#ifndef Eikonal_H_INC
#define Eikonal_H_INC



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
 * @return  Returns the number of cubes in the blob
 * @param[in/out] Distance  Distance function
 * @param[in] ID            Segmentation id
 * @param[in] DM            Domain information
 * @param[in] timesteps     Maximum number of timesteps to process
 * @return  Returns the global variation
 */
inline float Eikonal3D( Array<float> &Distance, const Array<char> &ID, const Domain &Dm, const int timesteps);


/*!
 * @brief  Calculate the distance using a simple method
 * @details  This routine calculates the distance using a very simple method working off the segmentation id.
 *
 * @return  Returns the number of cubes in the blob
 * @param[in/out] Distance  Distance function
 * @param[in] ID            Segmentation id
 * @param[in] DM            Domain information
 * @return  Returns the global variation
 */
inline void CalcDist3D( Array<float> &Distance, const Array<char> &ID, const Domain &Dm );


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
inline void CalcDistMultiLevel( Array<float> &Distance, const Array<char> &ID, const Domain &Dm );



#include "analysis/eikonal.hpp"

#endif

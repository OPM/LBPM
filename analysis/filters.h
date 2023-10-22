#ifndef Filters_H_INC
#define Filters_H_INC

#include "common/Array.h"

/*!
 * @brief  Filter image
 * @details  This routine performs a mean filter
 * @param[in] Input     Input image
 * @param[out] Output   Output image
 */
void Mean3D(const Array<double> &Input, Array<double> &Output);

/*!
 * @brief  Filter image
 * @details  This routine performs a median filter
 * @param[in] Input     Input image
 * @param[out] Output   Output image
 */
void Med3D(const Array<float> &Input, Array<float> &Output);

/*!
 * @brief  Filter image
 * @details  This routine performs a non-linear local means filter
 * @param[in] Input     Input image
 * @param[in] Mean      Mean value
 * @param[in] Distance  Distance
 * @param[out] Output   Output image
 * @param[in] d
 * @param[in] h
 */
int NLM3D(const Array<float> &Input, Array<float> &Mean,
          const Array<float> &Distance, Array<float> &Output, const int d,
          const float h);

#endif

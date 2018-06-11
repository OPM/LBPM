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
#ifndef Filters_H_INC
#define Filters_H_INC


#include "common/Array.h"


/*!
 * @brief  Filter image
 * @details  This routine performs a median filter
 * @param[in] Input     Input image
 * @param[out] Output   Output image
 */
void Med3D( const Array<float> &Input, Array<float> &Output );


/*!
 * @brief  Filter image
 * @details  This routine performs a non-linear local means filter
 * @param[in] Input     Input image
 * @param[in] Mean      Mean value
 * @param[out] Output   Output image
 */
int NLM3D( const Array<float> &Input, Array<float> &Mean, 
    const Array<float> &Distance, Array<float> &Output, const int d, const float h);


#endif

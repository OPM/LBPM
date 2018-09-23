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
#ifndef Distance_H_INC
#define Distance_H_INC

#include "common/Domain.h"
#include "common/Array.hpp"


struct Vec {
    double x;
    double y;
    double z;
    inline Vec(): x(0), y(0), z(0) {}
    inline Vec( double x_, double y_, double z_ ): x(x_), y(y_), z(z_) {}
    inline double norm() const { return sqrt(x*x+y*y+z*z); }
    inline double norm2() const { return x*x+y*y+z*z; }
};
inline bool operator<(const Vec& l, const Vec& r){ return l.x*l.x+l.y*l.y+l.z*l.z < r.x*r.x+r.y*r.y+r.z*r.z; }


/*!
 * @brief  Calculate the distance using a simple method
 * @details  This routine calculates the vector distance to the nearest domain surface.
 * @param[out] Distance     Distance function
 * @param[in] ID            Segmentation id
 * @param[in] Dm            Domain information
 * @param[in] periodic      Directions that are periodic
 */
template<class TYPE>
void CalcDist( Array<TYPE> &Distance, const Array<char> &ID, const Domain &Dm,
    const std::array<bool,3>& periodic = {true,true,true}, const std::array<double,3>& dx = {1,1,1} );

/*!
 * @brief  Calculate the distance using a simple method
 * @details  This routine calculates the vector distance to the nearest domain surface.
 * @param[out] Distance     Distance function
 * @param[in] ID            Domain id
 * @param[in] Dm            Domain information
 * @param[in] periodic      Directions that are periodic
 */
void CalcVecDist( Array<Vec> &Distance, const Array<int> &ID, const Domain &Dm,
    const std::array<bool,3>& periodic = {true,true,true}, const std::array<double,3>& dx = {1,1,1} );

#endif

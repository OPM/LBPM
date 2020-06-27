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

inline double minmod(double &a, double &b){

  double value;

  value = a;
  if ( a*b < 0.0)    value=0.0;
  else if (fabs(a) > fabs(b)) value = b;

  return value;
}

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


/*!
 * @brief  Calculate the distance based on solution of Eikonal equation
 * @details  This routine calculates the signed distance to the nearest domain surface.
 * @param[out] Distance     Distance function
 * @param[in] ID            Domain id
 * @param[in] Dm            Domain information
 * @param[in] timesteps      number of timesteps to run for Eikonal solver
 * @param[in] periodic      Directions that are periodic 
 */
double Eikonal(DoubleArray &Distance, char *ID, Domain &Dm, int timesteps, const std::array<bool,3>& periodic);

#endif

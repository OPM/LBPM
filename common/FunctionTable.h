/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University
  Copyright Equnior ASA

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
#ifndef included_FunctionTable
#define included_FunctionTable

#include "common/ArraySize.h"

#include <functional>

/*!
 * Class FunctionTable is a serial function table class that defines
 *   a series of operations that can be performed on the Array class.
 *   Users can impliment additional versions of the function table that match
 *   the interface to change the behavior of the array class.
 */
class FunctionTable final {
public:
    /*!
     * Initialize the array with random values
     * @param[in] x         The array to operate on
     */
    template <class TYPE, class FUN> static void rand(Array<TYPE, FUN> &x);

    /*!
     * Perform a reduce operator y = f(x)
     * @param[in] op            The function operation
     *                          Note: the operator is a template parameter to improve performance
     * @param[in] A             The array to operate on
     * @param[in] initialValue  The initial value for the reduction (0 for sum, +/- inf for min/max,
     * ...)
     * @return                  The reduction
     */
    template <class TYPE, class FUN, typename LAMBDA>
    static inline TYPE reduce(LAMBDA &op, const Array<TYPE, FUN> &A,
                              const TYPE &initialValue);

    /*!
     * Perform a reduce operator z = f(x,y)
     * @param[in] op            The function operation
     *                          Note: the operator is a template parameter to improve performance
     * @param[in] A             The first array to operate on
     * @param[in] B             The second array to operate on
     * @param[in] initialValue  The initial value for the reduction (0 for sum, +/- inf for min/max,
     * ...)
     * @return                  The reduction
     */
    template <class TYPE, class FUN, typename LAMBDA>
    static inline TYPE reduce(LAMBDA &op, const Array<TYPE, FUN> &A,
                              const Array<TYPE, FUN> &B,
                              const TYPE &initialValue);

    /*!
     * Perform a element-wise operation y = f(x)
     * @param[in] fun           The function operation
     *                          Note: the function is a template parameter to improve performance
     * @param[in,out] x         The array to operate on
     * @param[out] y            The output array
     */
    template <class TYPE, class FUN, typename LAMBDA>
    static inline void transform(LAMBDA &fun, const Array<TYPE, FUN> &x,
                                 Array<TYPE, FUN> &y);

    /*!
     * Perform a element-wise operation z = f(x,y)
     * @param[in] fun           The function operation
     *                          Note: the function is a template parameter to improve performance
     * @param[in] x             The first array
     * @param[in] y             The second array
     * @param[out] z            The output array
     */
    template <class TYPE, class FUN, typename LAMBDA>
    static inline void transform(LAMBDA &fun, const Array<TYPE, FUN> &x,
                                 const Array<TYPE, FUN> &y,
                                 Array<TYPE, FUN> &z);

    /*!
     * Multiply two arrays
     * @param[in] a             The first array
     * @param[in] b             The second array
     * @param[out] c            The output array
     */
    template <class TYPE, class FUN>
    static void multiply(const Array<TYPE, FUN> &a, const Array<TYPE, FUN> &b,
                         Array<TYPE, FUN> &c);

    /*!
     * Perform dgemv/dgemm equavalent operation ( C = alpha*A*B + beta*C )
     * @param[in] alpha         The scalar value alpha
     * @param[in] A             The first array
     * @param[in] B             The second array
     * @param[in] beta          The scalar value alpha
     * @param[in,out] C         The output array C
     */
    template <class TYPE, class FUN>
    static void gemm(const TYPE alpha, const Array<TYPE, FUN> &A,
                     const Array<TYPE, FUN> &B, const TYPE beta,
                     Array<TYPE, FUN> &C);

    /*!
     * Perform axpy equavalent operation ( y = alpha*x + y )
     * @param[in] alpha         The scalar value alpha
     * @param[in] x             The input array x
     * @param[in,out] y         The output array y
     */
    template <class TYPE, class FUN>
    static void axpy(const TYPE alpha, const Array<TYPE, FUN> &x,
                     Array<TYPE, FUN> &y);

    /*!
     * Check if two arrays are approximately equal
     * @param[in] A             The first array
     * @param[in] B             The second array
     * @param[in] tol           The tolerance
     */
    template <class TYPE, class FUN>
    static bool equals(const Array<TYPE, FUN> &A, const Array<TYPE, FUN> &B,
                       TYPE tol);

    template <class TYPE>
    static inline void gemmWrapper(char TRANSA, char TRANSB, int M, int N,
                                   int K, TYPE alpha, const TYPE *A, int LDA,
                                   const TYPE *B, int LDB, TYPE beta, TYPE *C,
                                   int LDC);

    /* Specialized Functions */

    /*!
     * Perform a element-wise operation y = max(x , 0)
     * @param[in] A             The input array
     * @param[out] B            The output array
     */
    template <class TYPE, class FUN, class ALLOC>
    static void transformReLU(const Array<TYPE, FUN, ALLOC> &A,
                              Array<TYPE, FUN, ALLOC> &B);

    /*!
     * Perform a element-wise operation B = |A|
     * @param[in] A             The array to operate on
     * @param[out] B            The output array
     */
    template <class TYPE, class FUN, class ALLOC>
    static void transformAbs(const Array<TYPE, FUN, ALLOC> &A,
                             Array<TYPE, FUN, ALLOC> &B);

    /*!
     * Perform a element-wise operation B = tanh(A)
     * @param[in] A             The array to operate on
     * @param[out] B            The output array
     */
    template <class TYPE, class FUN, class ALLOC>
    static void transformTanh(const Array<TYPE, FUN, ALLOC> &A,
                              Array<TYPE, FUN, ALLOC> &B);

    /*!
     * Perform a element-wise operation B = max(-1 , min(1 , A) )
     * @param[in] A             The array to operate on
     * @param[out] B            The output array
     */
    template <class TYPE, class FUN, class ALLOC>
    static void transformHardTanh(const Array<TYPE, FUN, ALLOC> &A,
                                  Array<TYPE, FUN, ALLOC> &B);

    /*!
     * Perform a element-wise operation B = 1 / (1 + exp(-A))
     * @param[in] A             The array to operate on
     * @param[out] B            The output array
     */
    template <class TYPE, class FUN, class ALLOC>
    static void transformSigmoid(const Array<TYPE, FUN, ALLOC> &A,
                                 Array<TYPE, FUN, ALLOC> &B);

    /*!
     * Perform a element-wise operation B = log(exp(A) + 1)
     * @param[in] A             The array to operate on
     * @param[out] B            The output array
     */
    template <class TYPE, class FUN, class ALLOC>
    static void transformSoftPlus(const Array<TYPE, FUN, ALLOC> &A,
                                  Array<TYPE, FUN, ALLOC> &B);

    /*!
     * Sum the elements of the Array
     * @param[in] A             The array to sum
     */
    template <class TYPE, class FUN, class ALLOC>
    static TYPE sum(const Array<TYPE, FUN, ALLOC> &A);

private:
    FunctionTable();

    template <class T> static inline void rand(size_t N, T *x);
};

#endif

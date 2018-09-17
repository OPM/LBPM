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
class FunctionTable final
{
public:
    /*!
     * Initialize the array with random values
     * @param[in] x         The array to operate on
     */
    template<class TYPE, class FUN>
    static void rand( Array<TYPE, FUN> &x );

    /*!
     * Perform a reduce operator y = f(x)
     * @param[in] op        The function operation
     *                      Note: the operator is a template parameter
     *                      (compared to a std::function to improve performance)
     * @param[in] A         The array to operate on
     * @return              The reduction
     */
    template<class TYPE, class FUN, typename LAMBDA>
    static inline TYPE reduce( LAMBDA &op, const Array<TYPE, FUN> &A );

    /*!
     * Perform a element-wise operation y = f(x)
     * @param[in] fun       The function operation
     *                      Note: the operator is a template parameter
     *                      (compared to a std::function to improve performance)
     * @param[in] x         The input array to operate on
     * @param[out] y        The output array
     */
    template<class TYPE, class FUN, typename LAMBDA>
    static inline void transform( LAMBDA &fun, const Array<TYPE, FUN> &x, Array<TYPE, FUN> &y );

    /*!
     * Perform a element-wise operation z = f(x,y)
     * @param[in] fun       The function operation
     *                      Note: the operator is a template parameter
     *                      (compared to a std::function to improve performance)
     * @param[in] x         The first array
     * @param[in] y         The second array
     * @param[out] z        The result
     */
    template<class TYPE, class FUN, typename LAMBDA>
    static inline void transform(
        LAMBDA &fun, const Array<TYPE, FUN> &x, const Array<TYPE, FUN> &y, Array<TYPE, FUN> &z );

    /*!
     * Multiply two arrays
     * @param[in] a             The first array
     * @param[in] b             The second array
     * @param[out] c            The output array
     */
    template<class TYPE, class FUN>
    static inline void multiply(
        const Array<TYPE, FUN> &a, const Array<TYPE, FUN> &b, Array<TYPE, FUN> &c );

    /*!
     * Perform dgemv/dgemm equavalent operation ( C = alpha*A*B + beta*C )
     * @param[in] alpha         The scalar value alpha
     * @param[in] A             The first array
     * @param[in] B             The second array
     * @param[in] beta          The scalar value alpha
     * @param[in,out] c         The output array C
     */
    template<class TYPE, class FUN>
    static void gemm( const TYPE alpha, const Array<TYPE, FUN> &A, const Array<TYPE, FUN> &B,
        const TYPE beta, Array<TYPE, FUN> &C );

    /*!
     * Perform axpy equavalent operation ( y = alpha*x + y )
     * @param[in] alpha         The scalar value alpha
     * @param[in] x             The input array x
     * @param[in,out] y         The output array y
     */
    template<class TYPE, class FUN>
    static void axpy( const TYPE alpha, const Array<TYPE, FUN> &x, Array<TYPE, FUN> &y );

    /*!
     * Check if two arrays are approximately equal
     * @param[in] A             The first array
     * @param[in] B             The second array
     * @param[in] tol           The tolerance
     */
    template<class TYPE, class FUN>
    static bool equals( const Array<TYPE, FUN> &A, const Array<TYPE, FUN> &B, TYPE tol );


private:
    FunctionTable();
};


#endif

#ifndef included_FunctionTable
#define included_FunctionTable


#include "common/Array.h"

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
    static void multiply(
        const Array<TYPE, FUN> &a, const Array<TYPE, FUN> &b, Array<TYPE, FUN> &c );


private:
    FunctionTable();

    template<class T>
    static inline void rand( size_t N, T *x );
};

#include "common/FunctionTable.hpp"

#endif

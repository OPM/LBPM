#ifndef included_ArraySizeClass
#define included_ArraySizeClass

#include "common/Utilities.h"

#include <array>
#include <cstdint>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <initializer_list>
#include <stdexcept>
#include <vector>

#if defined(__CUDA_ARCH__)
#include <cuda.h>
#define HOST_DEVICE __host__ __device__
#else
#define HOST_DEVICE
#endif
#if defined(USING_GCC) || defined(USING_CLANG)
#define ARRAY_ATTRIBUTE HOST_DEVICE __attribute__((always_inline))
#else
#define ARRAY_ATTRIBUTE HOST_DEVICE
#endif

#if (defined(DEBUG) || defined(_DEBUG)) && !defined(NDEBUG)
#define CHECK_ARRAY_LENGTH(i, length)                                          \
    do {                                                                       \
        if (i >= length)                                                       \
            throw std::out_of_range("Index exceeds array bounds");             \
    } while (0)
#else
#define CHECK_ARRAY_LENGTH(i, length)                                          \
    do {                                                                       \
    } while (0)
#endif

// Forward declerations
class FunctionTable;
template <class TYPE, class FUN = FunctionTable,
          class Allocator = std::allocator<TYPE>>
class Array;

//! Simple range class
template <class TYPE = size_t> class Range final {
public:
    //! Empty constructor
    Range() : i(0), j(-1), k(1) {}

    /*!
     * Create a range i:k:j (or i:j)
     * @param i_            Starting value
     * @param j_            Ending value
     * @param k_            Increment value
     */
    Range(const TYPE &i_, const TYPE &j_, const TYPE &k_ = 1)
        : i(i_), j(j_), k(k_) {}

    //! Get the number of values in the range
    size_t size() const {
        if (std::is_integral<TYPE>::value) {
            return (static_cast<int64_t>(j) - static_cast<int64_t>(i)) /
                   static_cast<int64_t>(k);
        } else if (std::is_floating_point<TYPE>::value) {
            double tmp = static_cast<double>((j - i)) / static_cast<double>(k);
            return static_cast<size_t>(floor(tmp + 1e-12) + 1);
        } else if (std::is_same<TYPE, std::complex<float>>::value ||
                   std::is_same<TYPE, std::complex<double>>::value) {
            double tmp = std::real((j - i) / (k));
            return static_cast<size_t>(floor(tmp + 1e-12) + 1);
        } else {
            ERROR("Unsupported type for range");
        }
    }

public:
    TYPE i, j, k;
};

//! Simple class to store the array dimensions
class ArraySize final {
public:
    //! Empty constructor
    ArraySize() : d_ndim(1), d_length(0), d_N{0, 1, 1, 1, 1} {}

    /*!
     * Create the vector size
     * @param N1            Number of elements in the first dimension
     */
    ArraySize(size_t N1) : d_ndim(1), d_length(N1), d_N{N1, 1, 1, 1, 1} {}

    /*!
     * Create the vector size
     * @param N1            Number of elements in the first dimension
     * @param N2            Number of elements in the second dimension
     */
    ArraySize(size_t N1, size_t N2)
        : d_ndim(2), d_length(N1 * N2), d_N{N1, N2, 1, 1, 1} {}

    /*!
     * Create the vector size
     * @param N1            Number of elements in the first dimension
     * @param N2            Number of elements in the second dimension
     * @param N3            Number of elements in the third dimension
     */
    ArraySize(size_t N1, size_t N2, size_t N3)
        : d_ndim(3), d_length(N1 * N2 * N3), d_N{N1, N2, N3, 1, 1} {}

    /*!
     * Create the vector size
     * @param N1            Number of elements in the first dimension
     * @param N2            Number of elements in the second dimension
     * @param N3            Number of elements in the third dimension
     * @param N4            Number of elements in the fourth dimension
     */
    ArraySize(size_t N1, size_t N2, size_t N3, size_t N4)
        : d_ndim(4), d_length(N1 * N2 * N3 * N4), d_N{N1, N2, N3, N4, 1} {}

    /*!
     * Create the vector size
     * @param N1            Number of elements in the first dimension
     * @param N2            Number of elements in the second dimension
     * @param N3            Number of elements in the third dimension
     * @param N4            Number of elements in the fourth dimension
     * @param N5            Number of elements in the fifth dimension
     */
    ArraySize(size_t N1, size_t N2, size_t N3, size_t N4, size_t N5)
        : d_ndim(5), d_length(N1 * N2 * N3 * N4 * N5), d_N{N1, N2, N3, N4, N5} {
    }

    /*!
     * Create from initializer list
     * @param N             Size of the array
     * @param ndim          Number of dimensions
     */
    ArraySize(std::initializer_list<size_t> N, int ndim = -1)
        : d_ndim(N.size()), d_length(0), d_N{0, 1, 1, 1, 1} {
        if (ndim >= 0)
            d_ndim = ndim;
        if (d_ndim > 5)
            throw std::out_of_range("Maximum number of dimensions exceeded");
        auto it = N.begin();
        for (size_t i = 0; i < d_ndim; i++, ++it)
            d_N[i] = *it;
        d_length = 1;
        for (unsigned long i : d_N)
            d_length *= i;
        if (d_ndim == 0)
            d_length = 0;
    }

    /*!
     * Create from raw pointer
     * @param ndim          Number of dimensions
     * @param dims          Dimensions
     */
    ArraySize(size_t ndim, const size_t *dims)
        : d_ndim(ndim), d_length(0), d_N{0, 1, 1, 1, 1} {
        if (d_ndim > 5)
            throw std::out_of_range("Maximum number of dimensions exceeded");
        for (size_t i = 0; i < ndim; i++)
            d_N[i] = dims[i];
        d_length = 1;
        for (unsigned long i : d_N)
            d_length *= i;
        if (d_ndim == 0)
            d_length = 0;
    }

    /*!
     * Create from std::array
     * @param N             Size of the array
     */
    template <std::size_t NDIM>
    ArraySize(const std::array<size_t, NDIM> &N) : ArraySize(NDIM, N.data()) {}

    /*!
     * Create from std::vector
     * @param N             Size of the array
     */
    inline ArraySize(const std::vector<size_t> &N)
        : ArraySize(N.size(), N.data()) {}

    // Copy/assignment constructors
    ArraySize(ArraySize &&rhs) = default;
    ArraySize(const ArraySize &rhs) = default;
    ArraySize &operator=(ArraySize &&rhs) = default;
    ArraySize &operator=(const ArraySize &rhs) = default;

    /*!
     * Access the ith dimension
     * @param i             Index to access
     */
    ARRAY_ATTRIBUTE size_t operator[](size_t i) const { return d_N[i]; }

    //! Return the number of dimensions
    ARRAY_ATTRIBUTE uint8_t ndim() const { return d_ndim; }

    //! Return the number of dimensions
    ARRAY_ATTRIBUTE size_t size() const { return d_ndim; }

    //! Return the total number of elements in the array
    ARRAY_ATTRIBUTE size_t length() const { return d_length; }

    //! Resize the dimension
    void resize(uint8_t dim, size_t N) {
        if (dim >= d_ndim)
            throw std::out_of_range("Invalid dimension");
        d_N[dim] = N;
        d_length = 1;
        for (unsigned long i : d_N)
            d_length *= i;
    }

    /*!
     * Reshape the Array so that the number of dimensions is the
     *    max of ndim and the largest dim>1.
     * @param ndim          Desired number of dimensions
     */
    void setNdim(uint8_t ndim) { d_ndim = std::max(ndim, d_ndim); }

    /*!
     * Remove singleton dimensions
     */
    void squeeze() {
        d_ndim = 0;
        for (uint8_t i = 0; i < maxDim(); i++) {
            if (d_N[i] != 1)
                d_N[d_ndim++] = d_N[i];
        }
    }

    //! Returns an iterator to the beginning
    const size_t *begin() const { return d_N; }

    //! Returns an iterator to the end
    const size_t *end() const { return d_N + d_ndim; }

    // Check if two array sizes are equal
    ARRAY_ATTRIBUTE bool operator==(const ArraySize &rhs) const {
        return d_ndim == rhs.d_ndim && memcmp(d_N, rhs.d_N, sizeof(d_N)) == 0;
    }

    // Check if two array sizes are equal (ignoring the dimension)
    ARRAY_ATTRIBUTE bool approxEqual(const ArraySize &rhs) const {
        return (length() == 0 && rhs.length() == 0) ||
               memcmp(d_N, rhs.d_N, sizeof(d_N)) == 0;
    }

    //! Check if two matrices are not equal
    ARRAY_ATTRIBUTE bool operator!=(const ArraySize &rhs) const {
        return d_ndim != rhs.d_ndim || memcmp(d_N, rhs.d_N, sizeof(d_N)) != 0;
    }

    //! Maximum supported dimension
    ARRAY_ATTRIBUTE static uint8_t maxDim() { return 5; }

    //! Get the index
    ARRAY_ATTRIBUTE size_t index(size_t i) const {
        CHECK_ARRAY_LENGTH(i, d_length);
        return i;
    }

    //! Get the index
    ARRAY_ATTRIBUTE size_t index(size_t i1, size_t i2) const {
        size_t index = i1 + i2 * d_N[0];
        CHECK_ARRAY_LENGTH(index, d_length);
        return index;
    }

    //! Get the index
    ARRAY_ATTRIBUTE size_t index(size_t i1, size_t i2, size_t i3) const {
        size_t index = i1 + d_N[0] * (i2 + d_N[1] * i3);
        CHECK_ARRAY_LENGTH(index, d_length);
        return index;
    }

    //! Get the index
    ARRAY_ATTRIBUTE size_t index(size_t i1, size_t i2, size_t i3,
                                 size_t i4) const {
        size_t index = i1 + d_N[0] * (i2 + d_N[1] * (i3 + d_N[2] * i4));
        CHECK_ARRAY_LENGTH(index, d_length);
        return index;
    }

    //! Get the index
    ARRAY_ATTRIBUTE size_t index(size_t i1, size_t i2, size_t i3, size_t i4,
                                 size_t i5) const {
        size_t index =
            i1 + d_N[0] * (i2 + d_N[1] * (i3 + d_N[2] * (i4 + d_N[3] * i5)));
        CHECK_ARRAY_LENGTH(index, d_length);
        return index;
    }

    //! Get the index
    size_t index(const std::array<size_t, 5> &i) const {
        size_t j = 0;
        for (size_t m = 0, N = 1; m < 5; m++) {
            j += i[m] * N;
            N *= d_N[m];
        }
        return j;
    }

    //! Get the index
    size_t index(std::initializer_list<size_t> i) const {
        size_t N = 1;
        size_t j = 0;
        size_t m = 0;
        for (size_t k : i) {
            j += k * N;
            N *= d_N[m++];
        }
        return j;
    }

    //! Convert the index to ijk values
    std::array<size_t, 5> ijk(size_t index) const {
        CHECK_ARRAY_LENGTH(index, d_length);
        size_t i0 = index % d_N[0];
        index = index / d_N[0];
        size_t i1 = index % d_N[1];
        index = index / d_N[1];
        size_t i2 = index % d_N[2];
        index = index / d_N[2];
        size_t i3 = index % d_N[3];
        index = index / d_N[3];
        return {i0, i1, i2, i3, index};
    }

    //! Convert the index to ijk values
    void ijk(size_t index, size_t *x) const {
        CHECK_ARRAY_LENGTH(index, d_length);
        x[0] = index % d_N[0];
        index = index / d_N[0];
        x[1] = index % d_N[1];
        index = index / d_N[1];
        x[2] = index % d_N[2];
        index = index / d_N[2];
        x[3] = index % d_N[3];
        index = index / d_N[3];
        x[4] = index;
    }

private:
    uint8_t d_ndim;
    size_t d_length;
    size_t d_N[5];
};

// Function to concatenate dimensions of two array sizes
inline ArraySize cat(const ArraySize &x, const ArraySize &y) {
    if (x.ndim() + y.ndim() > 5)
        throw std::out_of_range("Maximum number of dimensions exceeded");
    size_t N[5] = {0};
    for (int i = 0; i < x.ndim(); i++)
        N[i] = x[i];
    for (int i = 0; i < y.ndim(); i++)
        N[i + x.ndim()] = y[i];
    return ArraySize(x.ndim() + y.ndim(), N);
}

// Operator overloads
inline ArraySize operator*(size_t v, const ArraySize &x) {
    size_t N[5] = {v * x[0], v * x[1], v * x[2], v * x[3], v * x[4]};
    return ArraySize(x.ndim(), N);
}
inline ArraySize operator*(const ArraySize &x, size_t v) {
    size_t N[5] = {v * x[0], v * x[1], v * x[2], v * x[3], v * x[4]};
    return ArraySize(x.ndim(), N);
}
inline ArraySize operator-(const ArraySize &x, size_t v) {
    size_t N[5] = {x[0] - v, x[1] - v, x[2] - v, x[3] - v, x[4] - v};
    return ArraySize(x.ndim(), N);
}
inline ArraySize operator+(const ArraySize &x, size_t v) {
    size_t N[5] = {x[0] + v, x[1] + v, x[2] + v, x[3] + v, x[4] + v};
    return ArraySize(x.ndim(), N);
}
inline ArraySize operator+(size_t v, const ArraySize &x) {
    size_t N[5] = {x[0] + v, x[1] + v, x[2] + v, x[3] + v, x[4] + v};
    return ArraySize(x.ndim(), N);
}

#if defined(USING_ICC)
ENABLE_WARNINGS
#endif

#endif

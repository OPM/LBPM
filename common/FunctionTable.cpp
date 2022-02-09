#include "FunctionTable.hpp"

/********************************************************
 *  Random number generation                             *
 ********************************************************/
/*template<> char genRand<char>()
{
    static std::random_device rd;
    static std::mt19937 gen( rd() );
    static std::uniform_int_distribution<char> dis;
    return dis( gen );
}
template<> int8_t genRand<int8_t>()
{
    static std::random_device rd;
    static std::mt19937 gen( rd() );
    static std::uniform_int_distribution<int8_t> dis;
    return dis( gen );
}
template<> uint8_t genRand<uint8_t>()
{
    static std::random_device rd;
    static std::mt19937 gen( rd() );
    static std::uniform_int_distribution<uint8_t> dis;
    return dis( gen );
}
template<> int16_t genRand<int16_t>()
{
    static std::random_device rd;
    static std::mt19937 gen( rd() );
    static std::uniform_int_distribution<int16_t> dis;
    return dis( gen );
}
template<> uint16_t genRand<uint16_t>()
{
    static std::random_device rd;
    static std::mt19937 gen( rd() );
    static std::uniform_int_distribution<uint16_t> dis;
    return dis( gen );
}
template<> int32_t genRand<int32_t>()
{
    static std::random_device rd;
    static std::mt19937 gen( rd() );
    static std::uniform_int_distribution<int32_t> dis;
    return dis( gen );
}
template<> uint32_t genRand<uint32_t>()
{
    static std::random_device rd;
    static std::mt19937 gen( rd() );
    static std::uniform_int_distribution<uint32_t> dis;
    return dis( gen );
}
template<> int64_t genRand<int64_t>()
{
    static std::random_device rd;
    static std::mt19937 gen( rd() );
    static std::uniform_int_distribution<int64_t> dis;
    return dis( gen );
}
template<> uint64_t genRand<uint64_t>()
{
    static std::random_device rd;
    static std::mt19937 gen( rd() );
    static std::uniform_int_distribution<uint64_t> dis;
    return dis( gen );
}
template<> float genRand<float>()
{
    static std::random_device rd;
    static std::mt19937 gen( rd() );
    static std::uniform_real_distribution<float> dis;
    return dis( gen );
}
template<> double genRand<double>()
{
    static std::random_device rd;
    static std::mt19937 gen( rd() );
    static std::uniform_real_distribution<double> dis;
    return dis( gen );
}
template<> long double genRand<long double>()
{
    static std::random_device rd;
    static std::mt19937 gen( rd() );
    static std::uniform_real_distribution<double> dis;
    return dis( gen );
}
*/

/********************************************************
 *  axpy                                                 *
 ********************************************************/
template <>
void call_axpy<float>(size_t, const float, const float*, float*) {
    ERROR("Not finished");
}
template <>
void call_axpy<double>(size_t, const double, const double*, double*) {
    ERROR("Not finished");
}

/********************************************************
 *  Multiply two arrays                                  *
 ********************************************************/
template <>
void call_gemv<double>(size_t, size_t, double, double,
                       const double*, const double*, double*) {
    ERROR("Not finished");
}
template <>
void call_gemv<float>(size_t, size_t, float, float,
                      const float*, const float*, float*) {
    ERROR("Not finished");
}
template <>
void call_gemm<double>(size_t, size_t, size_t, double, double,
                       const double*, const double*, double*) {
    ERROR("Not finished");
}
template <>
void call_gemm<float>(size_t, size_t, size_t, float, float,
                      const float*, const float*, float*) {
    ERROR("Not finished");
}

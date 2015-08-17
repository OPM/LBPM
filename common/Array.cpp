#include "Array.h"
#include <stdint.h>


/********************************************************
*  std::swap                                            *
********************************************************/
namespace std
{
    template<> void swap( Array<bool>& v1, Array<bool>& v2 ) { v1.swap(v2); }
    template<> void swap( Array<char>& v1, Array<char>& v2 ) { v1.swap(v2); }
    template<> void swap( Array<int>& v1, Array<int>& v2 ) { v1.swap(v2); }
    template<> void swap( Array<unsigned int>& v1, Array<unsigned int>& v2 ) { v1.swap(v2); }
    template<> void swap( Array<int64_t>& v1, Array<int64_t>& v2 ) { v1.swap(v2); }
    template<> void swap( Array<uint64_t>& v1, Array<uint64_t>& v2 ) { v1.swap(v2); }
    template<> void swap( Array<float>& v1, Array<float>& v2 ) { v1.swap(v2); }
    template<> void swap( Array<double>& v1, Array<double>& v2 ) { v1.swap(v2); }
}



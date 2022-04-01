// This file contains unctions to pack/unpack data structures
#ifndef included_PackData
#define included_PackData

#include <map>
#include <set>
#include <vector>
#include <cstddef>

//! Template function to return the buffer size required to pack a class
template<class TYPE>
size_t packsize( const TYPE &rhs );

//! Template function to pack a class to a buffer
template<class TYPE>
void pack( const TYPE &rhs, char *buffer );

//! Template function to unpack a class from a buffer
template<class TYPE>
void unpack( TYPE &data, const char *buffer );


//! Template function to return the buffer size required to pack a std::vector
template<class TYPE>
size_t packsize( const std::vector<TYPE> &rhs );

//! Template function to pack a class to a buffer
template<class TYPE>
void pack( const std::vector<TYPE> &rhs, char *buffer );

//! Template function to pack a class to a buffer
template<class TYPE>
void unpack( std::vector<TYPE> &data, const char *buffer );


//! Template function to return the buffer size required to pack a std::pair
template<class TYPE1, class TYPE2>
size_t packsize( const std::pair<TYPE1, TYPE2> &rhs );

//! Template function to pack a class to a buffer
template<class TYPE1, class TYPE2>
void pack( const std::pair<TYPE1, TYPE2> &rhs, char *buffer );

//! Template function to pack a class to a buffer
template<class TYPE1, class TYPE2>
void unpack( std::pair<TYPE1, TYPE2> &data, const char *buffer );


//! Template function to return the buffer size required to pack a std::map
template<class TYPE1, class TYPE2>
size_t packsize( const std::map<TYPE1, TYPE2> &rhs );

//! Template function to pack a class to a buffer
template<class TYPE1, class TYPE2>
void pack( const std::map<TYPE1, TYPE2> &rhs, char *buffer );

//! Template function to pack a class to a buffer
template<class TYPE1, class TYPE2>
void unpack( std::map<TYPE1, TYPE2> &data, const char *buffer );


//! Template function to return the buffer size required to pack a std::set
template<class TYPE>
size_t packsize( const std::set<TYPE> &rhs );

//! Template function to pack a class to a buffer
template<class TYPE>
void pack( const std::set<TYPE> &rhs, char *buffer );

//! Template function to pack a class to a buffer
template<class TYPE>
void unpack( std::set<TYPE> &data, const char *buffer );


#include "IO/PackData.hpp"

#endif

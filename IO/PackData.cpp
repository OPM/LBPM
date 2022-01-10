#include "IO/PackData.h"
#include <string>


/********************************************************
 * Concrete implimentations for packing/unpacking        *
 ********************************************************/
// unsigned char
template<>
size_t packsize<unsigned char>( const unsigned char &rhs )
{
    return sizeof( unsigned char );
}
template<>
void pack<unsigned char>( const unsigned char &rhs, char *buffer )
{
    memcpy( buffer, &rhs, sizeof( unsigned char ) );
}
template<>
void unpack<unsigned char>( unsigned char &data, const char *buffer )
{
    memcpy( &data, buffer, sizeof( unsigned char ) );
}
// char
template<>
size_t packsize<char>( const char &rhs )
{
    return sizeof( char );
}
template<>
void pack<char>( const char &rhs, char *buffer )
{
    memcpy( buffer, &rhs, sizeof( char ) );
}
template<>
void unpack<char>( char &data, const char *buffer )
{
    memcpy( &data, buffer, sizeof( char ) );
}
// int
template<>
size_t packsize<int>( const int &rhs )
{
    return sizeof( int );
}
template<>
void pack<int>( const int &rhs, char *buffer )
{
    memcpy( buffer, &rhs, sizeof( int ) );
}
template<>
void unpack<int>( int &data, const char *buffer )
{
    memcpy( &data, buffer, sizeof( int ) );
}
// unsigned int
template<>
size_t packsize<unsigned int>( const unsigned int &rhs )
{
    return sizeof( unsigned int );
}
template<>
void pack<unsigned int>( const unsigned int &rhs, char *buffer )
{
    memcpy( buffer, &rhs, sizeof( int ) );
}
template<>
void unpack<unsigned int>( unsigned int &data, const char *buffer )
{
    memcpy( &data, buffer, sizeof( int ) );
}
// size_t
template<>
size_t packsize<size_t>( const size_t &rhs )
{
    return sizeof( size_t );
}
template<>
void pack<size_t>( const size_t &rhs, char *buffer )
{
    memcpy( buffer, &rhs, sizeof( size_t ) );
}
template<>
void unpack<size_t>( size_t &data, const char *buffer )
{
    memcpy( &data, buffer, sizeof( size_t ) );
}
// std::string
template<>
size_t packsize<std::string>( const std::string &rhs )
{
    return rhs.size() + 1;
}
template<>
void pack<std::string>( const std::string &rhs, char *buffer )
{
    memcpy( buffer, rhs.c_str(), rhs.size() + 1 );
}
template<>
void unpack<std::string>( std::string &data, const char *buffer )
{
    data = std::string( buffer );
}

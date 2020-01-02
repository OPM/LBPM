#ifndef included_StackTrace_stringView
#define included_StackTrace_stringView

#include <cstring>
#include <ostream>

namespace StackTrace {

// string_view
class string_view
{
public:
    // Constants:
    static constexpr size_t npos = size_t( -1 );

    // Constructions
    constexpr string_view() noexcept : d_data( nullptr ), d_size( 0 ) {}
    constexpr string_view( string_view&& ) noexcept      = default;
    constexpr string_view( const string_view& ) noexcept = default;
    constexpr string_view( const char* s ) : d_data( s ), d_size( s ? strlen( s ) : 0 ) {}
    constexpr string_view( const char* s, size_t count ) : d_data( s ), d_size( count ) {}
    inline string_view( const std::string& s ) : d_data( s.data() ), d_size( s.size() ) {}

    // Assignment
    constexpr string_view& operator=( string_view&& other ) noexcept = default;
    constexpr string_view& operator=( const string_view& other ) noexcept = default;

    // Iterators
    constexpr const char* begin() const noexcept { return d_data; }
    constexpr const char* end() const noexcept { return d_data + d_size; }
    constexpr const char* cbegin() const noexcept { return begin(); }
    constexpr const char* cend() const noexcept { return end(); }

    // capacity
    constexpr size_t size() const noexcept { return d_size; }
    constexpr size_t length() const noexcept { return d_size; }
    constexpr bool empty() const noexcept { return d_size == 0; }

    // Element access
    constexpr const char& operator[]( size_t pos ) const
    {
        if ( pos >= d_size )
            throw std::out_of_range( "string_view[]" );
        return d_data[pos];
    }
    constexpr const char& at( size_t pos ) const
    {
        if ( pos >= d_size )
            throw std::out_of_range( "string_view::at()" );
        return d_data[pos];
    }
    constexpr const char& front() const
    {
        if ( d_size == 0 )
            throw std::out_of_range( "front()" );
        return d_data[0];
    }
    constexpr const char& back() const
    {
        if ( d_size == 0 )
            throw std::out_of_range( "back()" );
        return d_data[size() - 1];
    }
    constexpr const char* data() const noexcept { return d_data; }

    // Swap data
    void swap( string_view& other ) noexcept
    {
        std::swap( d_data, other.d_data );
        std::swap( d_size, other.d_size );
    }

    // String operations
    size_t copy( char* dest, size_t n, size_t pos = 0 ) const
    {
        if ( pos > size() )
            throw std::out_of_range( "string_view::copy()" );
        const size_t rlen = std::min( n, size() - pos );
        memcpy( dest, data() + pos, rlen );
        return rlen;
    }
    constexpr string_view substr( size_t pos = 0, size_t n = npos ) const
    {
        if ( pos > size() )
            throw std::out_of_range( "string_view::substr()" );
        return string_view( data() + pos, std::min( n, size() - pos ) );
    }

    // Find
    constexpr size_t find( char ch, size_t pos = 0 ) const noexcept
    {
        for ( size_t i = pos; i < d_size; i++ )
            if ( d_data[i] == ch )
                return i;
        return std::string::npos;
    }
    constexpr size_t find( string_view v, size_t pos = 0 ) const noexcept
    {
        size_t i = pos;
        size_t N = v.size();
        if ( N == 0 || N > ( d_size - pos ) )
            return std::string::npos;
        while ( i < ( d_size - N + 1 ) ) {
            size_t j = 0;
            for ( j = 0; j < N && i + j < d_size; j++ )
                if ( d_data[i + j] != v[j] )
                    break;
            if ( j == N )
                return i;
            i++;
        }
        return std::string::npos;
    }

    // compare()
    constexpr int compare( const string_view& other ) const noexcept
    {
        int N      = std::min( size(), other.size() );
        int result = 0;
        for ( int i = 0; i < N && result == 0; i++ )
            if ( d_data[i] != other[i] )
                result = d_data[i] < other[i] ? -( i + 1 ) : ( i + 1 );
        if ( result == 0 )
            result = size() == other.size() ? 0 : size() < other.size() ? -1 : 1;
        return result;
    }
    constexpr int compare( size_t pos1, size_t n1, string_view other ) const
    {
        return substr( pos1, n1 ).compare( other );
    }
    constexpr int compare( size_t pos1, size_t n1, string_view other, size_t pos2, size_t n2 ) const
    {
        return substr( pos1, n1 ).compare( other.substr( pos2, n2 ) );
    }
    constexpr int compare( char const* s ) const { return compare( string_view( s ) ); }
    constexpr int compare( size_t pos1, size_t n1, char const* s ) const
    {
        return substr( pos1, n1 ).compare( string_view( s ) );
    }
    constexpr int compare( size_t pos1, size_t n1, char const* s, size_t n2 ) const
    {
        return substr( pos1, n1 ).compare( string_view( s, n2 ) );
    }

    explicit operator std::string() const { return std::string( begin(), end() ); }
    std::string to_string() const { return std::string( begin(), end() ); }

private:
    const char* d_data;
    size_t d_size;
};


// Non-member functions:
constexpr inline bool operator==( const string_view& lhs, const string_view& rhs ) noexcept
{
    return lhs.compare( rhs ) == 0;
}
constexpr inline bool operator!=( const string_view& lhs, const string_view& rhs ) noexcept
{
    return lhs.compare( rhs ) != 0;
}
constexpr inline bool operator<( const string_view& lhs, const string_view& rhs ) noexcept
{
    return lhs.compare( rhs ) < 0;
}

constexpr inline bool operator<=( const string_view& lhs, const string_view& rhs ) noexcept
{
    return lhs.compare( rhs ) <= 0;
}
constexpr inline bool operator>( const string_view& lhs, const string_view& rhs ) noexcept
{
    return lhs.compare( rhs ) > 0;
}
constexpr inline bool operator>=( const string_view& lhs, const string_view& rhs ) noexcept
{
    return lhs.compare( rhs ) >= 0;
}
inline std::string to_string( const string_view& v ) { return std::string( v.begin(), v.end() ); }
inline string_view to_string_view( std::string const& s )
{
    return string_view( s.data(), s.size() );
}
inline std::ostream& operator<<( std::ostream& out, const string_view& s )
{
    out << s.data();
    return out;
}

} // namespace StackTrace

#endif

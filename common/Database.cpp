/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University

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
#include "common/Database.h"
#include "common/Utilities.h"

#include <algorithm>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <string>
#include <tuple>


/********************************************************************
 * Constructors/destructor                                           *
 ********************************************************************/
Database::Database()  = default;
Database::~Database() = default;
Database::Database( const Database& rhs ) : KeyData( rhs )
{
    d_data.clear();
    for ( const auto& tmp : rhs.d_data )
        putData( tmp.first, tmp.second->clone() );
}
Database& Database::operator=( const Database& rhs )
{
    if ( this == &rhs )
        return *this;
    d_data.clear();
    for ( const auto& tmp : rhs.d_data )
        putData( tmp.first, tmp.second->clone() );
    return *this;
}
Database::Database( Database&& rhs ) { std::swap( d_data, rhs.d_data ); }
Database& Database::operator=( Database&& rhs )
{
    if ( this != &rhs )
        std::swap( d_data, rhs.d_data );
    return *this;
}


/********************************************************************
 * Clone the database                                                *
 ********************************************************************/
std::shared_ptr<KeyData> Database::clone() const { return cloneDatabase(); }
std::shared_ptr<Database> Database::cloneDatabase() const
{
    auto db = std::make_shared<Database>();
    for ( const auto& tmp : d_data )
        db->putData( tmp.first, tmp.second->clone() );
    return db;
}


/********************************************************************
 * Get the data object                                               *
 ********************************************************************/
bool Database::keyExists( const std::string& key ) const
{
    return d_data.find( key ) != d_data.end();
}
std::shared_ptr<KeyData> Database::getData( const std::string& key )
{
    auto it = d_data.find( key );
    if ( it == d_data.end() ) {
        char msg[1000];
        sprintf( msg, "Variable %s was not found in database", key.c_str() );
        ERROR( msg );
    }
    return it->second;
}
std::shared_ptr<const KeyData> Database::getData( const std::string& key ) const
{
    return const_cast<Database*>( this )->getData( key );
}
bool Database::isDatabase( const std::string& key ) const
{
    auto ptr  = getData( key );
    auto ptr2 = std::dynamic_pointer_cast<const Database>( ptr );
    return ptr2 != nullptr;
}
std::shared_ptr<Database> Database::getDatabase( const std::string& key )
{
    std::shared_ptr<KeyData> ptr   = getData( key );
    std::shared_ptr<Database> ptr2 = std::dynamic_pointer_cast<Database>( ptr );
    if ( ptr2 == nullptr ) {
        char msg[1000];
        sprintf( msg, "Variable %s is not a database", key.c_str() );
        ERROR( msg );
    }
    return ptr2;
}
std::shared_ptr<const Database> Database::getDatabase( const std::string& key ) const
{
    return const_cast<Database*>( this )->getDatabase( key );
}
std::vector<std::string> Database::getAllKeys() const
{
    std::vector<std::string> keys;
    keys.reserve( d_data.size() );
    for ( const auto& it : d_data )
        keys.push_back( it.first );
    return keys;
}
void Database::putDatabase( const std::string& key, std::shared_ptr<Database> db )
{
    d_data[key] = std::move( db );
}
void Database::putData( const std::string& key, std::shared_ptr<KeyData> data )
{
    d_data[key] = std::move( data );
}


/********************************************************************
 * Is the data of the given type                                     *
 ********************************************************************/
template<>
bool Database::isType<double>( const std::string& key ) const
{
    auto type = getData( key )->type();
    return type == "double";
}
template<>
bool Database::isType<float>( const std::string& key ) const
{
    auto type = getData( key )->type();
    return type == "double";
}
template<>
bool Database::isType<int>( const std::string& key ) const
{
    bool pass = true;
    auto type = getData( key )->type();
    if ( type == "double" ) {
        auto data = getVector<double>( key );
        for ( auto tmp : data )
            pass = pass && static_cast<double>( static_cast<int>( tmp ) ) == tmp;
    } else {
        pass = false;
    }
    return pass;
}
template<>
bool Database::isType<std::string>( const std::string& key ) const
{
    auto type = getData( key )->type();
    return type == "string";
}
template<>
bool Database::isType<bool>( const std::string& key ) const
{
    auto type = getData( key )->type();
    return type == "bool";
}


/********************************************************************
 * Get a vector                                                      *
 ********************************************************************/
template<>
std::vector<std::string> Database::getVector<std::string>(
    const std::string& key, const Units& ) const
{
    std::shared_ptr<const KeyData> ptr = getData( key );
    if ( std::dynamic_pointer_cast<const EmptyKeyData>( ptr ) )
        return std::vector<std::string>();
    const auto* ptr2 = dynamic_cast<const KeyDataString*>( ptr.get() );
    if ( ptr2 == nullptr ) {
        ERROR( "Key '" + key + "' is not a string" );
    }
    return ptr2->d_data;
}
template<>
std::vector<bool> Database::getVector<bool>( const std::string& key, const Units& ) const
{
    std::shared_ptr<const KeyData> ptr = getData( key );
    if ( std::dynamic_pointer_cast<const EmptyKeyData>( ptr ) )
        return std::vector<bool>();
    const auto* ptr2 = dynamic_cast<const KeyDataBool*>( ptr.get() );
    if ( ptr2 == nullptr ) {
        ERROR( "Key '" + key + "' is not a bool" );
    }
    return ptr2->d_data;
}
template<class TYPE>
std::vector<TYPE> Database::getVector( const std::string& key, const Units& unit ) const
{
    std::shared_ptr<const KeyData> ptr = getData( key );
    if ( std::dynamic_pointer_cast<const EmptyKeyData>( ptr ) )
        return std::vector<TYPE>();
    std::vector<TYPE> data;
    if ( std::dynamic_pointer_cast<const KeyDataDouble>( ptr ) ) {
        const auto* ptr2                 = dynamic_cast<const KeyDataDouble*>( ptr.get() );
        const std::vector<double>& data2 = ptr2->d_data;
        double factor                    = 1;
        if ( !unit.isNull() ) {
            INSIST( !ptr2->d_unit.isNull(), "Field " + key + " must have units" );
            factor = ptr2->d_unit.convert( unit );
            INSIST( factor != 0, "Unit conversion failed" );
        }
        data.resize( data2.size() );
        for ( size_t i = 0; i < data2.size(); i++ )
            data[i] = static_cast<TYPE>( factor * data2[i] );
    } else if ( std::dynamic_pointer_cast<const KeyDataString>( ptr ) ) {
        ERROR( "Converting std::string to another type" );
    } else if ( std::dynamic_pointer_cast<const KeyDataBool>( ptr ) ) {
        ERROR( "Converting std::bool to another type" );
    } else {
        ERROR( "Unable to convert data format" );
    }
    return data;
}


/********************************************************************
 * Put a vector                                                      *
 ********************************************************************/
template<>
void Database::putVector<std::string>(
    const std::string& key, const std::vector<std::string>& data, const Units& )
{
    std::shared_ptr<KeyDataString> ptr( new KeyDataString() );
    ptr->d_data = data;
    d_data[key] = ptr;
}
template<>
void Database::putVector<bool>(
    const std::string& key, const std::vector<bool>& data, const Units& )
{
    std::shared_ptr<KeyDataBool> ptr( new KeyDataBool() );
    ptr->d_data = data;
    d_data[key] = ptr;
}
template<class TYPE>
void Database::putVector( const std::string& key, const std::vector<TYPE>& data, const Units& unit )
{
    std::shared_ptr<KeyDataDouble> ptr( new KeyDataDouble() );
    ptr->d_unit = unit;
    ptr->d_data.resize( data.size() );
    for ( size_t i = 0; i < data.size(); i++ )
        ptr->d_data[i] = static_cast<double>( data[i] );
    d_data[key] = ptr;
}


/********************************************************************
 * Print the database                                                *
 ********************************************************************/
void Database::print( std::ostream& os, const std::string& indent ) const
{
    for ( const auto& it : d_data ) {
        os << indent << it.first;
        if ( dynamic_cast<const Database*>( it.second.get() ) ) {
            const auto* db = dynamic_cast<const Database*>( it.second.get() );
            os << " {\n";
            db->print( os, indent + "   " );
            os << indent << "}\n";
        } else {
            os << " = ";
            it.second->print( os, "" );
        }
    }
}
std::string Database::print( const std::string& indent ) const
{
    std::stringstream ss;
    print( ss, indent );
    return ss.str();
}


/********************************************************************
 * Read input database file                                          *
 ********************************************************************/
Database::Database( const std::string& filename )
{
    // Read the input file into memory
    FILE* fid = fopen( filename.c_str(), "rb" );
    if ( fid == nullptr )
        ERROR( "Error opening file " + filename );
    fseek( fid, 0, SEEK_END );
    size_t bytes = ftell( fid );
    rewind( fid );
    auto* buffer  = new char[bytes + 4];
    size_t result = fread( buffer, 1, bytes, fid );
    fclose( fid );
    if ( result != bytes )
        ERROR( "Error reading file " + filename );
    buffer[bytes + 0] = '\n';
    buffer[bytes + 1] = '}';
    buffer[bytes + 2] = '\n';
    buffer[bytes + 3] = 0;
    // Create the database entries
    loadDatabase( buffer, *this );
    // Free temporary memory
    delete[] buffer;
}
std::shared_ptr<Database> Database::createFromString( const std::string& data )
{
    std::shared_ptr<Database> db( new Database() );
    auto* buffer = new char[data.size() + 4];
    memcpy( buffer, data.data(), data.size() );
    buffer[data.size() + 0] = '\n';
    buffer[data.size() + 1] = '}';
    buffer[data.size() + 2] = '\n';
    buffer[data.size() + 3] = 0;
    loadDatabase( buffer, *db );
    delete[] buffer;
    return db;
}
enum class token_type {
    newline,
    line_comment,
    block_start,
    block_stop,
    quote,
    equal,
    bracket,
    end_bracket,
    end
};
inline size_t length( token_type type )
{
    size_t len = 0;
    if ( type == token_type::newline || type == token_type::quote || type == token_type::equal ||
         type == token_type::bracket || type == token_type::end_bracket ||
         type == token_type::end ) {
        len = 1;
    } else if ( type == token_type::line_comment || type == token_type::block_start ||
                type == token_type::block_stop ) {
        len = 2;
    }
    return len;
}
inline std::tuple<size_t, token_type> find_next_token( const char* buffer )
{
    size_t i = 0;
    while ( true ) {
        if ( buffer[i] == '\n' || buffer[i] == '\r' ) {
            return std::pair<size_t, token_type>( i + 1, token_type::newline );
        } else if ( buffer[i] == 0 ) {
            return std::pair<size_t, token_type>( i + 1, token_type::end );
        } else if ( buffer[i] == '"' ) {
            return std::pair<size_t, token_type>( i + 1, token_type::quote );
        } else if ( buffer[i] == '=' ) {
            return std::pair<size_t, token_type>( i + 1, token_type::equal );
        } else if ( buffer[i] == '{' ) {
            return std::pair<size_t, token_type>( i + 1, token_type::bracket );
        } else if ( buffer[i] == '}' ) {
            return std::pair<size_t, token_type>( i + 1, token_type::end_bracket );
        } else if ( buffer[i] == '/' ) {
            if ( buffer[i + 1] == '/' ) {
                return std::pair<size_t, token_type>( i + 2, token_type::line_comment );
            } else if ( buffer[i + 1] == '*' ) {
                return std::pair<size_t, token_type>( i + 2, token_type::block_start );
            }
        } else if ( buffer[i] == '*' ) {
            if ( buffer[i + 1] == '/' )
                return std::pair<size_t, token_type>( i + 2, token_type::block_stop );
        }
        i++;
    }
    return std::pair<size_t, token_type>( 0, token_type::end );
}
inline std::string deblank( const std::string& str )
{
    size_t i1 = 0xFFFFFFF, i2 = 0;
    for ( size_t i = 0; i < str.size(); i++ ) {
        if ( str[i] != ' ' ) {
            i1 = std::min( i1, i );
            i2 = std::max( i2, i );
        }
    }
    return i1 <= i2 ? str.substr( i1, i2 - i1 + 1 ) : std::string();
}
size_t skip_comment( const char* buffer )
{
    auto tmp                     = find_next_token( buffer );
    const token_type end_comment = ( std::get<1>( tmp ) == token_type::line_comment ) ?
                                       token_type::newline :
                                       token_type::block_stop;
    size_t pos = 0;
    while ( std::get<1>( tmp ) != end_comment ) {
        if ( std::get<1>( tmp ) == token_type::end )
            ERROR( "Encountered end of file before block comment end" );
        pos += std::get<0>( tmp );
        tmp = find_next_token( &buffer[pos] );
    }
    pos += std::get<0>( tmp );
    return pos;
}
inline std::string lower( const std::string& str )
{
    std::string tmp( str );
    std::transform( tmp.begin(), tmp.end(), tmp.begin(), ::tolower );
    return tmp;
}
static std::tuple<size_t, std::shared_ptr<KeyData>> read_value(
    const char* buffer, const std::string& key )
{
    // Get the value as a std::string
    size_t pos            = 0;
    token_type type       = token_type::end;
    std::tie( pos, type ) = find_next_token( &buffer[pos] );
    size_t len            = pos - length( type );
    while ( type != token_type::newline ) {
        if ( type == token_type::quote ) {
            size_t i            = 0;
            std::tie( i, type ) = find_next_token( &buffer[pos] );
            pos += i;
            while ( type != token_type::quote ) {
                ASSERT( type != token_type::end );
                std::tie( i, type ) = find_next_token( &buffer[pos] );
                pos += i;
            }
        } else if ( type == token_type::line_comment || type == token_type::block_start ) {
            len = pos - length( type );
            pos += skip_comment( &buffer[pos - length( type )] ) - length( type );
            break;
        }
        size_t i            = 0;
        std::tie( i, type ) = find_next_token( &buffer[pos] );
        pos += i;
        len = pos - length( type );
    }
    const std::string value = deblank( std::string( buffer, len ) );
    // Split the value to an array of values
    std::vector<std::string> values;
    size_t i0 = 0, i = 0, count = 0;
    for ( ; i < value.size(); i++ ) {
        if ( value[i] == '"' ) {
            count++;
        } else if ( value[i] == ',' && count % 2 == 0 ) {
            values.push_back( deblank( value.substr( i0, i - i0 ) ) );
            i0 = i + 1;
        }
    }
    values.push_back( deblank( value.substr( i0 ) ) );
    // Convert the string value to the database value
    std::shared_ptr<KeyData> data;
    if ( value.empty() ) {
        data.reset( new EmptyKeyData() );
    } else if ( value.find( '"' ) != std::string::npos ) {
        auto* data2 = new KeyDataString();
        data.reset( data2 );
        data2->d_data.resize( values.size() );
        for ( size_t i = 0; i < values.size(); i++ ) {
            ASSERT( values[i].size() >= 2 );
            ASSERT( values[i][0] == '"' && values[i][values[i].size() - 1] == '"' );
            data2->d_data[i] = values[i].substr( 1, values[i].size() - 2 );
        }
    } else if ( lower( value ) == "true" || lower( value ) == "false" ) {
        auto* data2 = new KeyDataBool();
        data.reset( data2 );
        data2->d_data.resize( values.size() );
        for ( size_t i = 0; i < values.size(); i++ ) {
            ASSERT( values[i].size() >= 2 );
            if ( lower( values[i] ) != "true" && lower( values[i] ) != "false" )
                ERROR( "Error converting " + key + " to logical array" );
            data2->d_data[i] = lower( values[i] ) == "true";
        }
    } else { // if ( value.find('.')!=std::string::npos || value.find('e')!=std::string::npos ) {
        auto* data2 = new KeyDataDouble();
        data.reset( data2 );
        data2->d_data.resize( values.size(), 0 );
        for ( size_t i = 0; i < values.size(); i++ ) {
            Units unit;
            std::tie( data2->d_data[i], unit ) = KeyDataDouble::read( values[i] );
            if ( !unit.isNull() )
                data2->d_unit = unit;
        }
        //} else {
        //    ERROR("Unable to determine data type: "+value);
    }
    return std::tuple<size_t, std::shared_ptr<KeyData>>( pos, data );
}
size_t Database::loadDatabase( const char* buffer, Database& db )
{
    size_t pos = 0;
    while ( true ) {
        size_t i;
        token_type type;
        std::tie( i, type ) = find_next_token( &buffer[pos] );
        const std::string key =
            deblank( std::string( &buffer[pos], std::max<int>( i - length( type ), 1 ) - 1 ) );
        if ( type == token_type::line_comment || type == token_type::block_start ) {
            // Comment
            INSIST( key.empty(), "Key should be empty: " + key );
            pos += skip_comment( &buffer[pos] );
        } else if ( type == token_type::newline ) {
            INSIST( key.empty(), "Key should be empty: " + key );
            pos += i;
        } else if ( type == token_type::equal ) {
            // Reading key/value pair
            ASSERT( !key.empty() );
            pos += i;
            std::shared_ptr<KeyData> data;
            std::tie( i, data ) = read_value( &buffer[pos], key );
            ASSERT( data.get() != nullptr );
            db.d_data[key] = data;
            pos += i;
        } else if ( type == token_type::bracket ) {
            // Read database
            ASSERT( !key.empty() );
            pos += i;
            std::shared_ptr<Database> database( new Database() );
            pos += loadDatabase( &buffer[pos], *database );
            db.d_data[key] = database;
        } else if ( type == token_type::end_bracket ) {
            // Finished with the database
            pos += i;
            break;
        } else {
            ERROR( "Error loading data" );
        }
    }
    return pos;
}


/********************************************************************
 * Data type helper functions                                        *
 ********************************************************************/
void KeyDataDouble::print( std::ostream& os, const std::string& indent ) const
{
    os << indent;
    for ( size_t i = 0; i < d_data.size(); i++ ) {
        if ( i > 0 )
            os << ", ";
        if ( d_data[i] != d_data[i] ) {
            os << "nan";
        } else if ( d_data[i] == std::numeric_limits<double>::infinity() ) {
            os << "inf";
        } else if ( d_data[i] == -std::numeric_limits<double>::infinity() ) {
            os << "-inf";
        } else {
            os << std::setprecision( 12 ) << d_data[i];
        }
    }
    if ( !d_unit.isNull() )
        os << " " << d_unit.str();
    os << std::endl;
}
std::tuple<double, Units> KeyDataDouble::read( const std::string& str )
{
    std::string tmp = deblank( str );
    size_t index    = tmp.find( " " );
    if ( index != std::string::npos ) {
        return std::make_tuple(
            readValue( tmp.substr( 0, index ) ), Units( tmp.substr( index + 1 ) ) );
    } else {
        return std::make_tuple( readValue( tmp ), Units() );
    }
}
double KeyDataDouble::readValue( const std::string& str )
{
    const std::string tmp = lower( str );
    double data           = 0;
    if ( tmp == "inf" || tmp == "infinity" ) {
        data = std::numeric_limits<double>::infinity();
    } else if ( tmp == "-inf" || tmp == "-infinity" ) {
        data = -std::numeric_limits<double>::infinity();
    } else if ( tmp == "nan" ) {
        data = std::numeric_limits<double>::quiet_NaN();
    } else if ( tmp.find( '/' ) != std::string::npos ) {
        ERROR( "Error reading value" );
    } else {
        char* pos = nullptr;
        data      = strtod( tmp.c_str(), &pos );
        if ( static_cast<size_t>( pos - tmp.c_str() ) == tmp.size() + 1 )
            ERROR( "Error reading value" );
    }
    return data;
}


/********************************************************************
 * Instantiations                                                    *
 ********************************************************************/
template std::vector<char> Database::getVector<char>( const std::string&, const Units& ) const;
template std::vector<int> Database::getVector<int>( const std::string&, const Units& ) const;
template std::vector<size_t> Database::getVector<size_t>( const std::string&, const Units& ) const;
template std::vector<float> Database::getVector<float>( const std::string&, const Units& ) const;
template std::vector<double> Database::getVector<double>( const std::string&, const Units& ) const;
template void Database::putVector<char>(
    const std::string&, const std::vector<char>&, const Units& );
template void Database::putVector<int>( const std::string&, const std::vector<int>&, const Units& );
template void Database::putVector<size_t>(
    const std::string&, const std::vector<size_t>&, const Units& );
template void Database::putVector<float>(
    const std::string&, const std::vector<float>&, const Units& );
template void Database::putVector<double>(
    const std::string&, const std::vector<double>&, const Units& );
template bool Database::isType<int>( const std::string& ) const;
template bool Database::isType<float>( const std::string& ) const;
template bool Database::isType<double>( const std::string& ) const;
template bool Database::isType<std::string>( const std::string& ) const;

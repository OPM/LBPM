#ifndef included_Database_hpp
#define included_Database_hpp

#include "common/Database.h"
#include "common/Utilities.h"

#include <tuple>


/********************************************************************
 * Basic classes for primative data types                            *
 ********************************************************************/
class EmptyKeyData : public KeyData
{
public:
    EmptyKeyData() {}
    virtual ~EmptyKeyData() {}
    virtual std::shared_ptr<KeyData> clone() const override
    {
        return std::make_shared<EmptyKeyData>();
    }
    virtual void print( std::ostream& os, const std::string& = "" ) const override
    {
        os << std::endl;
    }
    virtual std::string type() const override { return ""; }
};
class KeyDataDouble : public KeyData
{
public:
    KeyDataDouble() {}
    explicit KeyDataDouble( const std::vector<double>& data, const Units& unit )
        : d_data( data ), d_unit( unit )
    {
    }
    virtual ~KeyDataDouble() {}
    virtual std::shared_ptr<KeyData> clone() const override
    {
        return std::make_shared<KeyDataDouble>( d_data, d_unit );
    }
    virtual void print( std::ostream& os, const std::string& indent = "" ) const override;
    virtual std::string type() const override { return "double"; }

    static std::tuple<double, Units> read( const std::string& );
    static double readValue( const std::string& );

public:
    std::vector<double> d_data;
    Units d_unit;
};
class KeyDataBool : public KeyData
{
public:
    KeyDataBool() {}
    explicit KeyDataBool( const std::vector<bool>& data ) : d_data( data ) {}
    virtual ~KeyDataBool() {}
    virtual std::shared_ptr<KeyData> clone() const override
    {
        return std::make_shared<KeyDataBool>( d_data );
    }
    virtual void print( std::ostream& os, const std::string& indent = "" ) const override
    {
        os << indent;
        for ( size_t i = 0; i < d_data.size(); i++ ) {
            if ( i > 0 ) {
                os << ", ";
            }
            if ( d_data[i] ) {
                os << "true";
            } else {
                os << "false";
            }
        }
        os << std::endl;
    }
    virtual std::string type() const override { return "bool"; }
    std::vector<bool> d_data;
};
class KeyDataString : public KeyData
{
public:
    KeyDataString() {}
    explicit KeyDataString( const std::vector<std::string>& data ) : d_data( data ) {}
    virtual ~KeyDataString() {}
    virtual std::shared_ptr<KeyData> clone() const override
    {
        return std::make_shared<KeyDataString>( d_data );
    }
    virtual void print( std::ostream& os, const std::string& indent = "" ) const override
    {
        os << indent;
        for ( size_t i = 0; i < d_data.size(); i++ ) {
            if ( i > 0 ) {
                os << ", ";
            }
            os << '"' << d_data[i] << '"';
        }
        os << std::endl;
    }
    virtual std::string type() const override { return "string"; }
    std::vector<std::string> d_data;
};


/********************************************************************
 * Get a vector                                                      *
 ********************************************************************/
template<class TYPE>
inline std::vector<TYPE> Database::getVector(
    const std::string& key, const std::string& unit ) const
{
    return getVector<TYPE>( key, Units( unit ) );
}
template<class TYPE>
inline void Database::putVector(
    const std::string& key, const std::vector<TYPE>& data, const std::string& unit )
{
    putVector<TYPE>( key, data, Units( unit ) );
}


/********************************************************************
 * Get a scalar                                                      *
 ********************************************************************/
template<class TYPE>
inline TYPE Database::getScalar( const std::string& key, const Units& unit ) const
{
    const std::vector<TYPE>& data = getVector<TYPE>( key, unit );
    if ( data.size() != 1 ) {
        char msg[1000];
        sprintf( msg, "Variable %s is not a scalar", key.c_str() );
        ERROR( msg );
    }
    return data[0];
}
template<class TYPE>
inline TYPE Database::getWithDefault(
    const std::string& key, const TYPE& value, const Units& unit ) const
{
    if ( !keyExists( key ) )
        return value;
    return getScalar<TYPE>( key, unit );
}
template<class TYPE>
inline void Database::putScalar( const std::string& key, const TYPE& data, const Units& unit )
{
    putVector<TYPE>( key, std::vector<TYPE>( 1, data ), unit );
}
template<class TYPE>
inline TYPE Database::getScalar( const std::string& key, const std::string& unit ) const
{
    return getScalar<TYPE>( key, Units( unit ) );
}
template<class TYPE>
inline TYPE Database::getWithDefault(
    const std::string& key, const TYPE& value, const std::string& unit ) const
{
    return getWithDefault<TYPE>( key, value, Units( unit ) );
}
template<class TYPE>
inline void Database::putScalar( const std::string& key, const TYPE& data, const std::string& unit )
{
    putScalar<TYPE>( key, data, Units( unit ) );
}
template<class TYPE>
inline void putVector(
    const std::string& key, const std::vector<TYPE>& data, const std::string& unit )
{
    putVector<TYPE>( key, data, Units( unit ) );
}


#endif

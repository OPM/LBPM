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
#ifndef included_Database
#define included_Database

#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "common/Units.h"


inline bool exists( const std::string& filename )
{
     std::ifstream domain( filename );
	 return domain.good();
}


//! Base class to hold data of a given type
class KeyData
{
protected:
    //! Empty constructor
    KeyData() {}

public:
    //! Destructor
    virtual ~KeyData() {}
    //! Copy the data
    virtual std::shared_ptr<KeyData> clone() const = 0;
    //! Print the data to a stream
    virtual void print( std::ostream& os, const std::string& indent = "" ) const = 0;
    //! Return the native data type
    virtual std::string type() const = 0;

protected:
    KeyData( const KeyData& ) {}
    KeyData& operator=( const KeyData& );
};


//! Class to a database
class Database : public KeyData
{
public:
    //! Empty constructor
    Database();

    /**
     * Open an database file.
     * @param filename       Name of input file to open
     */
    explicit Database( const std::string& filename );

    /**
     * Create database from string
     * @param data       String containing the database data
     */
    static std::shared_ptr<Database> createFromString( const std::string& data );

    //! Copy constructor
    Database( const Database& );

    //! Assignment operator
    Database& operator=( const Database& );

    //! Move constructor
    Database( Database&& rhs );

    //! Move assignment operator
    Database& operator=( Database&& rhs );

    //! Destructor
    virtual ~Database();

    //! Copy the data
    virtual std::shared_ptr<KeyData> clone() const override;

    //! Copy the data
    std::shared_ptr<Database> cloneDatabase() const;


    /**
     * Return true if the specified key exists in the database and false
     *     otherwise.
     * @param[in] key           Key name to lookup.
     */
    bool keyExists( const std::string& key ) const;


    /**
     * Return all keys in the database.
     */
    std::vector<std::string> getAllKeys() const;


    //! Return the number of entries in the database
    size_t size() const { return d_data.size(); }


    /**
     * Get the scalar entry from the database with the specified key
     * name.  If the specified key does not exist in the database or
     * is not a scalar of the given type, then an error message is printed and
     * the program exits.
     *
     * @param[in] key           Key name in database.
     * @param[in] unit          Desired units
     */
    template<class TYPE>
    inline TYPE getScalar( const std::string& key, const Units& unit = Units() ) const;


    /// @copydoc Database::getScalar(const std::string&,const Units&) const
    template<class TYPE>
    inline TYPE getScalar( const std::string& key, const std::string& unit ) const;


    /**
     * Get the scalar entry from the database with the specified key
     * name.  If the specified key does not exist in the database the
     * the default value will be printed
     *
     * @param[in] key           Key name in database
     * @param[in] value         Default value
     * @param[in] unit          Desired units
     */
    template<class TYPE>
    inline TYPE getWithDefault(
        const std::string& key, const TYPE& value, const Units& unit = Units() ) const;


    /// @copydoc Database::getWithDefault(const std::string&,const TYPE&,const Units&) const
    template<class TYPE>
    inline TYPE getWithDefault(
        const std::string& key, const TYPE& value, const std::string& unit ) const;


    /**
     * Put the scalar entry into the database with the specified key name.  
     * @param key           Key name in database.
     * @param value         Value to store
     * @param unit          Desired units
     */
    template<class TYPE>
    inline void putScalar( const std::string& key, const TYPE& value, const Units& unit = Units() );


    /**
     * Put the scalar entry into the database with the specified key name.  
     * @param key           Key name in database.
     * @param value         Value to store
     * @param unit          Desired units
     */
    template<class TYPE>
    inline void putScalar( const std::string& key, const TYPE& value, const std::string& unit );


    /**
     * Get the vector entries from the database with the specified key
     * name.  If the specified key does not exist in the database or
     * is not of the given type, then an error message is printed and
     * the program exits.
     *
     * @param key           Key name in database.
     * @param unit          Desired units
     */
    template<class TYPE>
    std::vector<TYPE> getVector( const std::string& key, const Units& unit = Units() ) const;


    /// @copydoc Database::getVector(const std::string&,const Units&) const
    template<class TYPE>
    inline std::vector<TYPE> getVector( const std::string& key, const std::string& unit ) const;


    /**
     * Put the vector entries into the database with the specified key
     * name.  If the specified key does not exist in the database or
     * is not of the given type, then an error message is printed and
     * the program exits.
     *
     * @param key           Key name in database.
     * @param data          Data to store
     * @param unit          Desired units
     */
    template<class TYPE>
    void putVector(
        const std::string& key, const std::vector<TYPE>& data, const Units& unit = Units() );


    /// @copydoc Database::putVector(const std::string&,const std::vector<TYPE>&,const Units&)
    template<class TYPE>
    inline void putVector(
        const std::string& key, const std::vector<TYPE>& data, const std::string& unit );


    /**
     * Get the data for a key in the database.  If the specified key
     * does not exist in the database an error message is printed and
     * the program exits.
     *
     * @param key Key name in database.
     */
    std::shared_ptr<KeyData> getData( const std::string& key );

    /**
     * Get the data for a key in the database.  If the specified key
     * does not exist in the database an error message is printed and
     * the program exits.
     *
     * @param key Key name in database.
     */
    std::shared_ptr<const KeyData> getData( const std::string& key ) const;


    /**
     * Put the data for a key in the database.
     *
     * @param key       Key name in database.
     * @param data      Data to store
     */
    void putData( const std::string& key, std::shared_ptr<KeyData> data );


    // Check if the key is a database object
    bool isDatabase( const std::string& key ) const;


    // Check if the entry can be stored as the given type
    template<class TYPE>
    bool isType( const std::string& key ) const;


    /**
     * Get the database for a key in the database.  If the specified key
     * does not exist in the database an error message is printed and
     * the program exits.
     *
     * @param key Key name in database.
     */
    std::shared_ptr<Database> getDatabase( const std::string& key );

    /**
     * Get the database for a key in the database.  If the specified key
     * does not exist in the database an error message is printed and
     * the program exits.
     *
     * @param key Key name in database.
     */
    std::shared_ptr<const Database> getDatabase( const std::string& key ) const;


    /**
     * Get the database for a key in the database.  If the specified key
     * does not exist in the database an error message is printed and
     * the program exits.
     *
     * @param key       Key name in database.
     * @param db        Database to store
     */
    void putDatabase( const std::string& key, std::shared_ptr<Database> db );


    /**
     * Print the data to a stream
     * @param os        Output stream
     * @param indent    Indenting to use before each line
     */
    virtual void print( std::ostream& os, const std::string& indent = "" ) const override;


    //! Print the type
    virtual std::string type() const override { return "database"; }


    /**
     * Print the data to a string
     * @return          Output string
     */
    std::string print( const std::string& indent = "" ) const;


protected:
    std::map<std::string, std::shared_ptr<KeyData>> d_data;

    // Function to load a database from a buffer
    static size_t loadDatabase( const char* buffer, Database& db );
};


#include "common/Database.hpp"

#endif

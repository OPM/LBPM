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
#include "common/Units.h"
#include "common/Utilities.h"

#include <algorithm>
#include <cmath>
#include <string>


constexpr double Units::d_pow10[22];
constexpr char Units::d_prefixSymbol[];


/********************************************************************
 * Constructors                                                      *
 ********************************************************************/
Units::Units() : d_prefix( UnitPrefix::unknown ), d_unit( UnitValue::unknown ) {}
Units::Units( UnitPrefix p, UnitValue u ) : d_prefix( p ), d_unit( u ) {}
Units::Units( const std::string& unit )
    : d_prefix( UnitPrefix::unknown ), d_unit( UnitValue::unknown )
{
    // Parse the string to get it into a more friendly format
    auto tmp = unit;
    tmp.erase( std::remove( tmp.begin(), tmp.end(), ' ' ), tmp.end() );
    // Check if the character '-' is present indicating a seperation between the prefix and unit
    size_t index = tmp.find( '-' );
    if ( index != std::string::npos ) {
        d_prefix = getUnitPrefix( tmp.substr( 0, index ) );
        d_unit   = getUnitValue( tmp.substr( index + 1 ) );
    } else {
        if ( tmp.size() <= 1 ) {
            d_prefix = UnitPrefix::none;
            d_unit   = getUnitValue( tmp );
        } else if ( tmp.substr( 0, 2 ) == "da" ) {
            d_prefix = UnitPrefix::deca;
            d_unit   = getUnitValue( tmp.substr( 2 ) );
        } else {
            d_prefix = getUnitPrefix( tmp.substr( 0, 1 ) );
            d_unit   = getUnitValue( tmp.substr( 1 ) );
            if ( d_prefix == UnitPrefix::unknown || d_unit == UnitValue::unknown ) {
                d_prefix = UnitPrefix::none;
                d_unit   = getUnitValue( tmp );
            }
        }
    }
}


/********************************************************************
 * Get prefix                                                        *
 ********************************************************************/
Units::UnitPrefix Units::getUnitPrefix( const std::string& str ) noexcept
{
    Units::UnitPrefix value = UnitPrefix::unknown;
    if ( str.empty() ) {
        value = UnitPrefix::none;
    } else if ( str == "yotta" || str == "Y" ) {
        value = UnitPrefix::yotta;
    } else if ( str == "zetta" || str == "Z" ) {
        value = UnitPrefix::zetta;
    } else if ( str == "exa" || str == "E" ) {
        value = UnitPrefix::exa;
    } else if ( str == "peta" || str == "P" ) {
        value = UnitPrefix::peta;
    } else if ( str == "tera" || str == "T" ) {
        value = UnitPrefix::tera;
    } else if ( str == "giga" || str == "G" ) {
        value = UnitPrefix::giga;
    } else if ( str == "mega" || str == "M" ) {
        value = UnitPrefix::mega;
    } else if ( str == "kilo" || str == "k" ) {
        value = UnitPrefix::kilo;
    } else if ( str == "hecto" || str == "h" ) {
        value = UnitPrefix::hecto;
    } else if ( str == "deca" || str == "da" ) {
        value = UnitPrefix::deca;
    } else if ( str == "deci" || str == "d" ) {
        value = UnitPrefix::deci;
    } else if ( str == "centi" || str == "c" ) {
        value = UnitPrefix::centi;
    } else if ( str == "milli" || str == "m" ) {
        value = UnitPrefix::milli;
    } else if ( str == "micro" || str == "u" ) {
        value = UnitPrefix::micro;
    } else if ( str == "nano" || str == "n" ) {
        value = UnitPrefix::nano;
    } else if ( str == "pico" || str == "p" ) {
        value = UnitPrefix::pico;
    } else if ( str == "femto" || str == "f" ) {
        value = UnitPrefix::femto;
    } else if ( str == "atto" || str == "a" ) {
        value = UnitPrefix::atto;
    } else if ( str == "zepto" || str == "z" ) {
        value = UnitPrefix::zepto;
    } else if ( str == "yocto" || str == "y" ) {
        value = UnitPrefix::yocto;
    }
    return value;
}


/********************************************************************
 * Get unit value                                                    *
 ********************************************************************/
Units::UnitValue Units::getUnitValue( const std::string& str ) noexcept
{
    Units::UnitValue value = UnitValue::unknown;
    if ( str == "meter" || str == "m" ) {
        value = UnitValue::meter;
    } else if ( str == "gram" || str == "g" ) {
        value = UnitValue::gram;
    } else if ( str == "second" || str == "s" ) {
        value = UnitValue::second;
    } else if ( str == "ampere" || str == "A" ) {
        value = UnitValue::ampere;
    } else if ( str == "kelvin" || str == "K" ) {
        value = UnitValue::kelvin;
    } else if ( str == "joule" || str == "J" ) {
        value = UnitValue::joule;
    } else if ( str == "ergs" || str == "erg" ) {
        value = UnitValue::erg;
    } else if ( str == "degree" || str == "degrees" ) {
        value = UnitValue::degree;
    } else if ( str == "radian" || str == "radians" ) {
        value = UnitValue::radian;
    }
    return value;
}


/********************************************************************
 * Get unit type                                                     *
 ********************************************************************/
Units::UnitType Units::getUnitType( UnitValue u ) noexcept
{
    switch ( u ) {
    case UnitValue::meter:
        return UnitType::length;
    case UnitValue::gram:
        return UnitType::mass;
    case UnitValue::second:
        return UnitType::time;
    case UnitValue::ampere:
        return UnitType::current;
    case UnitValue::kelvin:
        return UnitType::temperature;
    case UnitValue::joule:
    case UnitValue::erg:
        return UnitType::energy;
    case UnitValue::degree:
    case UnitValue::radian:
        return UnitType::angle;
    default:
        return UnitType::unknown;
    }
}


/********************************************************************
 * Convert to another unit system                                    *
 ********************************************************************/
double Units::convert( const Units& rhs ) const noexcept
{
    if ( this->operator==( rhs ) )
        return 1;
    // Convert the prefix
    double cp = convert( d_prefix ) / convert( rhs.d_prefix );
    if ( d_unit == rhs.d_unit )
        return cp; // Only need to convert prefix
    // Convert the unit
    if ( getUnitType( d_unit ) != getUnitType( rhs.d_unit ) )
        return 0; // Invalid conversion
    double cu = 0;
    if ( d_unit == UnitValue::joule && rhs.d_unit == UnitValue::erg )
        cu = 1e7;
    else if ( d_unit == UnitValue::erg && rhs.d_unit == UnitValue::joule )
        cu = 1e-7;
    else if ( d_unit == UnitValue::degree && rhs.d_unit == UnitValue::radian )
        cu = 0.017453292519943;
    else if ( d_unit == UnitValue::radian && rhs.d_unit == UnitValue::degree )
        cu = 57.295779513082323;
    // Return the total conversion
    return cp * cu;
}


/********************************************************************
 * Write a string for the units                                      *
 ********************************************************************/
std::string Units::str() const
{
    ASSERT( !isNull() );
    return std::string( str( d_prefix ).data() ) + str( d_unit );
}
std::array<char, 3> Units::str( UnitPrefix p ) noexcept
{
    std::array<char, 3> str;
    str[0] = d_prefixSymbol[static_cast<int8_t>( p )];
    str[1] = 0;
    str[2] = 0;
    if ( p == UnitPrefix::deca )
        str[1] = 'a';
    return str;
}
std::string Units::str( UnitValue u )
{
    if ( u == UnitValue::meter ) {
        return "m";
    } else if ( u == UnitValue::gram ) {
        return "g";
    } else if ( u == UnitValue::second ) {
        return "s";
    } else if ( u == UnitValue::ampere ) {
        return "A";
    } else if ( u == UnitValue::kelvin ) {
        return "K";
    } else if ( u == UnitValue::joule ) {
        return "J";
    } else if ( u == UnitValue::erg ) {
        return "erg";
    } else if ( u == UnitValue::degree ) {
        return "degree";
    } else if ( u == UnitValue::radian ) {
        return "radian";
    }
    return "unknown";
}

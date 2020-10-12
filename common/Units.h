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
#ifndef included_Units
#define included_Units

#include <array>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>


//! Unit system class
class Units final
{
public:
    //! Enum to hold prefix
    enum class UnitPrefix : int8_t {
        yocto   = 0,
        zepto   = 1,
        atto    = 2,
        femto   = 3,
        pico    = 4,
        nano    = 5,
        micro   = 6,
        milli   = 7,
        centi   = 8,
        deci    = 9,
        none    = 10,
        deca    = 11,
        hecto   = 12,
        kilo    = 13,
        mega    = 14,
        giga    = 15,
        tera    = 16,
        peta    = 17,
        exa     = 18,
        zetta   = 19,
        yotta   = 20,
        unknown = 21
    };

    //! Enum to unit type
    enum class UnitType : uint8_t {
        length,
        mass,
        time,
        current,
        temperature,
        energy,
        angle,
        unknown
    };

    //! Enum to hold unit
    enum class UnitValue : uint8_t {
        meter,
        gram,
        second,
        ampere,
        kelvin,
        joule,
        erg,
        degree,
        radian,
        unknown
    };


public:
    //! Constructor
    Units();

    //! Constructor
    explicit Units( const std::string& unit );

    //! Constructor
    explicit Units( UnitPrefix, UnitValue );

    //! Get the prefix
    inline UnitPrefix getPrefix() const noexcept { return d_prefix; }

    //! Get the unit
    inline UnitValue getUnit() const noexcept { return d_unit; }

    //! Get the unit
    inline UnitType getUnitType() const noexcept { return getUnitType( d_unit ); }

    //! Get the unit
    static UnitType getUnitType( UnitValue ) noexcept;

    //! Get the prefix from a string
    static UnitPrefix getUnitPrefix( const std::string& ) noexcept;

    //! Get the unit value from a string
    static UnitValue getUnitValue( const std::string& ) noexcept;

    //! Convert to the given unit system
    double convert( const Units& ) const noexcept;

    //! Convert a prefix to a scalar
    static inline double convert( UnitPrefix x ) noexcept
    {
        return d_pow10[static_cast<int8_t>( x )];
    }

    //! Get a string representation of the units
    std::string str() const;

    //! Get a string representation for the prefix
    static std::array<char, 3> str( UnitPrefix ) noexcept;

    //! Get a string representation for the unit value
    static std::string str( UnitValue );

    //! Operator ==
    inline bool operator==( const Units& rhs ) const noexcept
    {
        return d_prefix == rhs.d_prefix && d_unit == rhs.d_unit;
    }

    //! Operator !=
    inline bool operator!=( const Units& rhs ) const noexcept
    {
        return d_prefix != rhs.d_prefix || d_unit != rhs.d_unit;
    }

    //! Check if unit is null
    bool isNull() const { return d_prefix == UnitPrefix::unknown || d_unit == UnitValue::unknown; }

protected:
    UnitPrefix d_prefix;
    UnitValue d_unit;

private:
    constexpr static double d_pow10[22]    = { 1e-24, 1e-21, 1e-18, 1e-15, 1e-12, 1e-9, 1e-6, 1e-3,
        1e-2, 0.1, 1, 10, 100, 1000, 1e6, 1e9, 1e12, 1e15, 1e18, 1e21, 1e24, 0 };
    constexpr static char d_prefixSymbol[] = "yzafpnumcd\0dhkMGTPEZYu";
};

#endif

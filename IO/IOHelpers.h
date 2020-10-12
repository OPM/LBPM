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
#ifndef IO_HELPERS_INC
#define IO_HELPERS_INC

#include <string.h>
#include <vector>

namespace IO {


// Find a character in a line
inline size_t find( const char *line, char token )
{
    size_t i=0;
    while ( 1 ) {
        if ( line[i]==token || line[i]<32 || line[i]==0 )
            break;
        ++i;
    }
    return i;
}


// Remove preceeding/trailing whitespace
inline std::string deblank( const std::string& str )
{
    size_t i1 = str.size();
    size_t i2 = 0;
    for (size_t i=0; i<str.size(); i++) {
        if ( str[i]!=' ' && str[i]>=32 ) {
            i1 = std::min(i1,i);
            i2 = std::max(i2,i);
        }
    }
    return str.substr(i1,i2-i1+1);
}


// Split a list by the given token
inline std::vector<std::string> splitList( const char *line, const char token )
{
    std::vector<std::string> list;
    size_t i1 = 0;
    size_t i2 = 0;
    while ( 1 ) {
        if ( line[i2]==token || line[i2]<32 ) {
            std::string tmp(&line[i1],i2-i1);
            tmp = deblank(tmp);
            if ( !tmp.empty() )
                list.push_back(tmp);
            i1 = i2+1;
        }
        if ( line[i2]==0 )
            break;
        i2++;
    }
    return list;
}



};

#endif


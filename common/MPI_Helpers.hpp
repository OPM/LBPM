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
// This file contains wrappers for MPI routines and functions to pack/unpack data structures
#ifndef MPI_WRAPPERS_HPP
#define MPI_WRAPPERS_HPP

#include "common/MPI_Helpers.h"
#include <string.h>
#include <vector>
#include <set>
#include <map>



/********************************************************
* Default instantiations for std::vector                *
********************************************************/
template<class TYPE>
size_t packsize( const std::vector<TYPE>& rhs )
{
    size_t bytes = sizeof(size_t);
    for (size_t i=0; i<rhs.size(); i++)
        bytes += packsize(rhs[i]);
    return bytes;
}
template<class TYPE>
void pack( const std::vector<TYPE>& rhs, char *buffer )
{
    size_t size = rhs.size();
    memcpy(buffer,&size,sizeof(size_t));
    size_t pos = sizeof(size_t);
    for (size_t i=0; i<rhs.size(); i++) {
        pack(rhs[i],&buffer[pos]);
        pos += packsize(rhs[i]);
    }
}
template<class TYPE>
void unpack( std::vector<TYPE>& data, const char *buffer )
{
    size_t size;
    memcpy(&size,buffer,sizeof(size_t));
    data.clear();
    data.resize(size);
    size_t pos = sizeof(size_t);
    for (size_t i=0; i<data.size(); i++) {
        unpack(data[i],&buffer[pos]);
        pos += packsize(data[i]);
    }
}


/********************************************************
* Default instantiations for std::pair                  *
********************************************************/
template<class TYPE1, class TYPE2>
size_t packsize( const std::pair<TYPE1,TYPE2>& rhs )
{
    return packsize(rhs.first)+packsize(rhs.second);
}
template<class TYPE1, class TYPE2>
void pack( const std::pair<TYPE1,TYPE2>& rhs, char *buffer )
{
    pack(rhs.first,buffer);
    pack(rhs.second,&buffer[packsize(rhs.first)]);
}
template<class TYPE1, class TYPE2>
void unpack( std::pair<TYPE1,TYPE2>& data, const char *buffer )
{
    unpack(data.first,buffer);
    unpack(data.second,&buffer[packsize(data.first)]);
}


/********************************************************
* Default instantiations for std::map                   *
********************************************************/
template<class TYPE1, class TYPE2>
size_t packsize( const std::map<TYPE1,TYPE2>& rhs )
{
    size_t bytes = sizeof(size_t);
    typename std::map<TYPE1,TYPE2>::const_iterator it;
    for (it=rhs.begin(); it!=rhs.end(); ++it) {
        bytes += packsize(it->first);
        bytes += packsize(it->second);
    }
    return bytes;
}
template<class TYPE1, class TYPE2>
void pack( const std::map<TYPE1,TYPE2>& rhs, char *buffer )
{
    size_t N = rhs.size();
    pack(N,buffer);
    size_t pos = sizeof(size_t);
    typename std::map<TYPE1,TYPE2>::const_iterator it;
    for (it=rhs.begin(); it!=rhs.end(); ++it) {
        pack(it->first,&buffer[pos]);   pos+=packsize(it->first);
        pack(it->second,&buffer[pos]);  pos+=packsize(it->second);
    }
}
template<class TYPE1, class TYPE2>
void unpack( std::map<TYPE1,TYPE2>& data, const char *buffer )
{
    size_t N = 0;
    unpack(N,buffer);
    size_t pos = sizeof(size_t);
    data.clear();
    for (size_t i=0; i<N; i++) {
        std::pair<TYPE1,TYPE2> tmp;
        unpack(tmp.first,&buffer[pos]);   pos+=packsize(tmp.first);
        unpack(tmp.second,&buffer[pos]);  pos+=packsize(tmp.second);
        data.insert(tmp);
    }
}


/********************************************************
* Default instantiations for std::set                   *
********************************************************/
template<class TYPE>
size_t packsize( const std::set<TYPE>& rhs )
{
    size_t bytes = sizeof(size_t);
    typename std::set<TYPE>::const_iterator it;
    for (it=rhs.begin(); it!=rhs.end(); ++it) {
        bytes += packsize(*it);
    }
    return bytes;
}
template<class TYPE>
void pack( const std::set<TYPE>& rhs, char *buffer )
{
    size_t N = rhs.size();
    pack(N,buffer);
    size_t pos = sizeof(size_t);
    typename std::set<TYPE>::const_iterator it;
    for (it=rhs.begin(); it!=rhs.end(); ++it) {
        pack(*it);   pos+=packsize(*it);
    }
}
template<class TYPE>
void unpack( std::set<TYPE>& data, const char *buffer )
{
    size_t N = 0;
    unpack(N,buffer);
    size_t pos = sizeof(size_t);
    data.clear();
    for (size_t i=0; i<N; i++) {
        TYPE tmp;
        unpack(tmp,&buffer[pos]);   pos+=packsize(tmp);
        data.insert(tmp);
    }
}


#endif


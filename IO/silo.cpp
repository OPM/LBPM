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
#include "IO/silo.h"
#include "common/Utilities.h"
#include "common/MPI_Helpers.h"

#include "ProfilerApp.h"


#ifdef USE_SILO

#include <silo.h>



namespace silo {


/****************************************************
* Open/close a file                                 *
****************************************************/
DBfile* open( const std::string& filename, FileMode mode )
{
    DBfile *fid = nullptr;
    if ( mode == CREATE ) {
        fid = DBCreate( filename.c_str(), DB_CLOBBER, DB_LOCAL, nullptr, DB_HDF5 );
    } else if ( mode == WRITE ) {
        fid = DBOpen( filename.c_str(), DB_HDF5, DB_APPEND );
    } else if ( mode == READ ) {
        fid = DBOpen( filename.c_str(), DB_HDF5, DB_READ );
    }
    return fid;
}
void close( DBfile* fid )
{
    DBClose( fid );
}


/****************************************************
* Helper functions                                  *
****************************************************/
VariableDataType varDataType( DBfile *fid, const std::string& name )
{
    auto type = DBGetVarType( fid, name.c_str() );
    VariableDataType type2 = VariableDataType::UNKNOWN;
    if ( type == DB_DOUBLE )
        type2 = VariableDataType::DOUBLE;
    else if ( type == DB_FLOAT )
        type2 = VariableDataType::FLOAT;
    else if ( type == DB_INT )
        type2 = VariableDataType::INT;
    return type2;
}


/****************************************************
* Write/read a uniform mesh to silo                 *
****************************************************/
void readUniformMesh( DBfile* fid, const std::string& meshname,
    std::vector<double>& range, std::vector<int>& N )
{
    DBquadmesh* mesh = DBGetQuadmesh( fid, meshname.c_str() );
    int ndim = mesh->ndims;
    range.resize(2*ndim);
    N.resize(ndim);
    for (int d=0; d<ndim; d++) {
        N[d] = mesh->dims[d]-1;
        range[2*d+0] = mesh->min_extents[d];
        range[2*d+1] = mesh->max_extents[d];
    }
    DBFreeQuadmesh( mesh );
}


/****************************************************
* Write a multimesh                                 *
****************************************************/
void writeMultiMesh( DBfile* fid, const std::string& meshname,
    const std::vector<std::string>& meshNames,
    const std::vector<int>& meshTypes )
{
    std::vector<char*> meshnames(meshNames.size());
    for ( size_t i = 0; i < meshNames.size(); ++i )
        meshnames[i] = (char *) meshNames[i].c_str();
    std::string tree_name = meshname + "_tree";
    DBoptlist *optList    = DBMakeOptlist( 1 );
    DBAddOption( optList, DBOPT_MRGTREE_NAME, (char *) tree_name.c_str() );
    DBPutMultimesh( fid, meshname.c_str(), meshNames.size(), meshnames.data(), (int*) meshTypes.data(), nullptr );
    DBFreeOptlist( optList );
}


/****************************************************
* Write a multivariable                             *
****************************************************/
void writeMultiVar( DBfile* fid, const std::string& varname,
    const std::vector<std::string>& varNames,
    const std::vector<int>& varTypes )
{
    std::vector<char*> varnames(varNames.size(),nullptr);
    for (size_t j=0; j<varNames.size(); j++)
        varnames[j] = const_cast<char*>(varNames[j].c_str());
    DBPutMultivar( fid, varname.c_str(), varNames.size(), varnames.data(), (int*) varTypes.data(), nullptr );
}



}; // silo namespace


#else

#endif

/*
 * Copyright (c) 2025 Diogo Nardelli Siebert <diogo.siebert@ufsc.br>
 *
 * Licensed under either of
 *   - Apache License, Version 2.0 (https://www.apache.org/licenses/LICENSE-2.0)
 *   - GNU General Public License, Version 3.0 or later (https://www.gnu.org/licenses/gpl-3.0.html)
 *
 * SPDX-License-Identifier: (Apache-2.0 OR GPL-3.0-or-later)
 */

#include "xmlvtk.h"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cstring>
#include <zlib.h>

using namespace std;

/** @brief Indicates whether data compression is used by default **/
bool Element::compress = false;
/** @brief Default cache size used during data processing. */
headerType cacheSize = 1500;

/**
 * @brief Formats a string using printf-style syntax.
 * @tparam Args Variadic arguments for formatting.
 * @param fmt Format string.
 * @return Formatted string.
 */
template<typename... Args>
std::string format_string(const char* fmt, Args&&... args)
{
    int size = std::snprintf(nullptr, 0, fmt, std::forward<Args>(args)...);
    if (size < 0) {
        throw std::runtime_error("format_string: snprintf error");
    }
    std::vector<char> buf(size + 1);
    int size2 = std::snprintf(buf.data(), buf.size(), fmt, std::forward<Args>(args)...);
    if (size2 < 0) {
        throw std::runtime_error("format_string: snprintf error");
    }
    
    return std::string(buf.data(), buf.data() + size2);
}

/**
 * @brief Encodes binary input to Base64 string representation.
 * @param input Input binary data.
 * @param len Length of input data.
 * @return Encoded Base64 string.
 */
std::string spc_base64_encode(const unsigned char* input, size_t len)
{
    static const char table[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    std::string output;
    output.reserve((len / 3 + (len % 3 != 0)) * 4);

    for (size_t i = 0; i < len; i += 3) {
        uint32_t n = (input[i] << 16) | (i + 1 < len ? input[i + 1] << 8 : 0) | (i + 2 < len ? input[i + 2] : 0);  // AAAAAABB BBBBCCCC CCDDDDDD
        output.push_back(table[(n >> 18) & 0x3F]);
        output.push_back(table[(n >> 12) & 0x3F]);
        output.push_back(i + 1 < len ? table[(n >> 6) & 0x3F] : '=');
        output.push_back(i + 2 < len ? table[n & 0x3F] : '=');
    }
    return output;
}

std::string CellData::header()
{
    return "<CellData>\n";
}

std::string CellData::footer()
{
    return "</CellData>\n";
}

std::string PointData::header()
{
    return "<PointData>\n";
}

std::string PointData::footer()
{
    return "</PointData>\n";
}

std::string AppendData::header()
{
    return "<AppendedData encoding=\"raw\">\n";
}

std::string AppendData::footer()
{
    return "</AppendedData>\n";
}

AppendData::AppendData()
{
    this -> totalSize = 0;
}

unsigned int AppendData::addData(unsigned char* pointer, unsigned int size)
{
    int offset = this->totalSize;
    this->sizeList.push_back(size);
    this->pointerList.push_back(pointer);
    this->totalSize += (size+4);
    return offset;
}

DataArray::DataArray(const std::string&  name_,const std::string&  type_, const std::string&  format_, int components_, unsigned char* pointer_, uint64_t points_)
	: components(components_), name(name_), type(type_), format(format_), pointer(pointer_), points(points_)

{
    if ( (this -> type == "Int8") || (this -> type == "UInt8") ) this -> typeSize = 1;
    else if ( (this -> type == "Int16") || (this -> type == "UInt16") ) this -> typeSize = 2;
    else if ( (this -> type == "Int32") || (this -> type == "UInt32") ||  (this -> type == "Float32") ) this -> typeSize = 4;
    else if ( (this -> type == "Int64") || (this -> type == "UInt64") ||  (this -> type == "Float64") ) this -> typeSize = 8;
    this -> dataSize = components * points;
}

std::string DataArray::header()
{
    std::ostringstream stringStream;
    stringStream << format_string("<DataArray Name=\"%s\" type=\"%s\" NumberOfComponents=\"%d\" format=\"%s\" %s", 
		    name.c_str(),  type.c_str() , components, format.c_str(), 
		    (format == "appended") ? format_string(" offset=\"%d\" />", offset).c_str() : ">");
    return stringStream.str();
}

std::string DataArray::footer()
{
    return "\n</DataArray>\n";
}

void VTIWriter::setCompress()
{
    this->compress = true;
}

PVTIWriter::PVTIWriter(const std::string& filename_)
	: filename(filename_), 
	  originX(0), originY(0), originZ(0),
	  wholeMinX(0), wholeMinY(0), wholeMinZ(0), wholeMaxX(0), wholeMaxY(0), wholeMaxZ(0), pieceCounter(0)
{
}

VTIWriter::VTIWriter(const std::string& filename_)
    : filename(filename_) ,  vtkVersion("0.1"),
      originX(0), originY(0), originZ(0)
{
    this -> compress = false;
    if (sizeof(headerType) == 8) this -> headerTypeName = "UInt64";
    else if (sizeof(headerType) == 4) this -> headerTypeName = "UInt32";
    this -> byteOrder = "LittleEndian";
}

void VTIWriter::addPointData(const std::string&  name,const std::string&  type, const std::string&  format, int components, unsigned char* pointer)
{
    std::unique_ptr<DataArray> data = std::make_unique<DataArray>( name , type , format , components , (unsigned char*) pointer ,  pd.points  );
    data -> compress = this -> compress;

    if (format == "appended")
    {
        data -> offset = ad.addData( data->pointer , data->dataSize * data->typeSize );
        ad.compress = this->compress;
    }

    pd.addChild( std::move(data) );
}

void VTIWriter::addCellData(const std::string&  name,const std::string&  type, const std::string&  format, int components, unsigned char* pointer)
{
    std::unique_ptr<DataArray> data = std::make_unique<DataArray>( name , type , format , components , (unsigned char*) pointer ,  cd.cells  );
    data -> compress = this -> compress;

    if (format == "appended")
    {
        data -> offset = ad.addData( data->pointer , data->dataSize * data->typeSize );
        ad.compress = this->compress;
    }

    cd.addChild( std::move(data) );
}

std::string VTIWriter::header()
{
    std::ostringstream stringStream;    
    stringStream << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << endl;
    stringStream << format_string("<VTKFile type=\"ImageData\" version=\"%s\" byte_order=\"%s\" header_type=\"%s\" %s>", "1.0", byteOrder.c_str() , headerTypeName.c_str() , compress ? "compressor=\"vtkZLibDataCompressor\"" : "" ) << endl;
    stringStream << format_string("<ImageData WholeExtent=\"%d %d %d %d %d %d\"", wholeMinX, wholeMaxX, wholeMinY, wholeMaxY, wholeMinZ, wholeMaxZ);
    stringStream << format_string(" Origin=\"%f %f %f\"", originX, originY, originZ);
    stringStream << format_string(" Spacing=\"%f %f %f \">", spaceX, spaceY, spaceZ) << endl;
    stringStream << format_string("<Piece Extent= \" %d %d %d %d %d %d\">" ,  pieceMinX , pieceMaxX , pieceMinY , pieceMaxY, pieceMinZ , pieceMaxZ) << endl;
    return stringStream.str();
}

std::string PVTIWriter::header()
{
    std::ostringstream stringStream;
    stringStream << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << endl;
    stringStream << "<VTKFile type=\"PImageData\" version=\"1.0\" >" <<  endl;
    stringStream << format_string("<PImageData WholeExtent=\"%d %d %d %d %d %d\"", wholeMinX, wholeMaxX, wholeMinY, wholeMaxY, wholeMinZ, wholeMaxZ);
    stringStream << format_string(" Origin=\"%f %f %f\"", originX, originY, originZ);
    stringStream << format_string(" Spacing=\"%f %f %f\">", spaceX, spaceY, spaceZ) << endl;
    
    return stringStream.str();
}

std::string VTIWriter::footer()
{
    return "</VTKFile>\n";
}

std::string PVTIWriter::footer()
{
    return "</PImageData>\n</VTKFile>\n";
}

void VTIWriter::setWholeExtent(int64_t minX_,int64_t minY_, int64_t minZ_,int64_t maxX_,int64_t maxY_, int64_t maxZ_)
{
    wholeMinX = minX_;
    wholeMinY = minY_;
    wholeMinZ = minZ_;
    wholeMaxX = maxX_;
    wholeMaxY = maxY_;
    wholeMaxZ = maxZ_;
    setPiece(minX_,minY_,minZ_,maxX_,maxY_, maxZ_);
}

std::ostream& operator<<(std::ostream& os, AppendData& obj)
{
    if (obj.totalSize > 0)
    {
        os << obj.header();
        os << "_";
        for (int k = 0; k < obj.pointerList.size(); k++)
        {
            if (obj.compress == false)
            {
                os.write( (char*) &obj.sizeList[k], sizeof(headerType) );
                os.write( (char*) obj.pointerList[k], obj.sizeList[k]  );
            }
            else
            {
                streampos beginPos = os.tellp();

                std::vector<unsigned char>compressedData (cacheSize);

                headerType totalByteSize = obj.sizeList[k];
                headerType numberOfBlocks =  totalByteSize / cacheSize  +  (totalByteSize % cacheSize != 0 );

                size_t infoSize = sizeof(headerType) * (3 + numberOfBlocks);
                std::vector<headerType>compressedInfo(infoSize );

                compressedInfo[0] = numberOfBlocks;
                compressedInfo[1] = (numberOfBlocks > 1) ?  cacheSize : totalByteSize;
                compressedInfo[2] = (totalByteSize % cacheSize == 0) ? cacheSize : totalByteSize % cacheSize ;

                os.write((char*) compressedInfo.data() , infoSize );

                for (int n = 1; n <= numberOfBlocks;  n++)
                {
                    uLongf compressedLength = cacheSize;
                    int numberOfBytesInBlock = (n < numberOfBlocks) ? compressedInfo[1] : compressedInfo[2];
                    const Bytef* src = reinterpret_cast<const Bytef*>(    obj.pointerList[k] + (n - 1) * cacheSize );
                    int ret = compress( (Bytef *) (compressedData.data() ), &compressedLength, src, numberOfBytesInBlock );
                    if (ret != Z_OK)
                    {
                        std::cerr << "Error compressing data , code =" << ret << std::endl;
                        throw std::runtime_error("Failed compressing data using zlib in AppendData");
                    }
                    compressedInfo[n+2] = compressedLength;

                    os.write( reinterpret_cast<const char*>(  compressedData.data() ) , compressedLength );
                }

                streampos endPos = os.tellp();

                os.seekp( beginPos );
                os.write( (char*) compressedInfo.data() , infoSize );
                os.seekp( endPos );
            }
        }        
        os << endl;
        os << obj.footer();

    }
    return os;
}

std::ostream& operator<<(std::ostream& os, DataArray& obj)
{
    os << obj.header() << endl;

    headerType totalByteSize = static_cast<headerType>( (obj.dataSize) * (obj.typeSize) );
    headerType count = min<headerType >( 12 - sizeof(headerType) , totalByteSize );

    if  (obj.format == "binary")
    {
        if (obj.compress == false)
        {
            os << spc_base64_encode( (unsigned char *) &totalByteSize, sizeof(headerType) );
            count = 0;

            int leftOverSize = 0;
            unsigned char leftOverBuffer[3];

            for (; count < totalByteSize ; count += cacheSize)
            {
                unsigned char* pointer = obj.pointer + count;
                size_t size = min<size_t>( cacheSize, totalByteSize - count);
                if (leftOverSize > 0)
                {
                        for (; leftOverSize < 3; leftOverSize++)
                        {
                            if (size > 0)
                            {
                                leftOverBuffer[leftOverSize] = *(pointer++);
                                size--;
                            }
                            else break;
                        }

                        os << spc_base64_encode( leftOverBuffer , leftOverSize );
                        leftOverSize = 0;
                }

                leftOverSize = (size % 3);
                size -= leftOverSize;
                for (int i = 0; i < leftOverSize; i++)
                {
                     leftOverBuffer[i] = pointer[size + i];
                }

                os << spc_base64_encode( pointer , size );

            }

            if (leftOverSize > 0)
            {
                os << spc_base64_encode( leftOverBuffer , leftOverSize );
                leftOverSize = 0;
            }
        }
        else
        {
            streampos beginPos = os.tellp();

            std::vector<unsigned char>compressedData (cacheSize + 3);
            int leftOver = 0;

            headerType numberOfBlocks =  totalByteSize / cacheSize  +  (totalByteSize % cacheSize != 0 );
            size_t infoSize = sizeof(headerType) * (3 + numberOfBlocks);
            std::vector<headerType>compressedInfo(infoSize );

            compressedInfo[0] = numberOfBlocks;
            compressedInfo[1] = (numberOfBlocks > 1) ?  cacheSize : totalByteSize;
            compressedInfo[2] = (totalByteSize % cacheSize == 0) ? cacheSize : totalByteSize % cacheSize ;

            os << spc_base64_encode( reinterpret_cast<unsigned char*>(compressedInfo.data()) , infoSize );

            for (int n = 1; n <= numberOfBlocks;  n++)
            {
                uLongf compressedLength = cacheSize;
                int numberOfBytesInBlock = (n < numberOfBlocks) ? compressedInfo[1] : compressedInfo[2];
                int ret = compress( reinterpret_cast<Bytef*>(compressedData.data() + 3) , &compressedLength, obj.pointer + (n-1) * cacheSize, numberOfBytesInBlock );
                if (ret != Z_OK)
                {
                    std::cerr << "Error compressing data , code =" << ret << std::endl;
                    throw std::runtime_error("Failed compressing data using zlib in AppendData");
                }
                compressedInfo[n+2] = compressedLength;
                
                int encodeSize = 3* ( (compressedLength + leftOver)/3 );
                os << spc_base64_encode( compressedData.data() + 3 - leftOver , encodeSize );

                int newLeftOver = compressedLength + leftOver - encodeSize ;
                if (newLeftOver > 0) 
                {
                    memcpy( compressedData.data() + 3 - newLeftOver, compressedData.data() + 3 - leftOver + encodeSize, newLeftOver);
                }
                leftOver = newLeftOver;
            }

            if (leftOver > 0)
            {
                os << spc_base64_encode(compressedData.data() + 3 - leftOver , leftOver );
            }

            streampos endPos = os.tellp();
            os.seekp( beginPos );

            os << spc_base64_encode( reinterpret_cast<unsigned char*>(compressedInfo.data()) , infoSize );
            os.seekp( endPos );
        }

        os << obj.footer();
    }
    
    return os;
}

void VTIWriter::setPiece(int64_t minX_,int64_t minY_, int64_t minZ_,int64_t maxX_,int64_t maxY_, int64_t maxZ_)
{
    pieceMinX = minX_;
    pieceMinY = minY_;
    pieceMinZ = minZ_;
    pieceMaxX = maxX_;
    pieceMaxY = maxY_;
    pieceMaxZ = maxZ_;
    sizeX = maxX_ - minX_ + 1;
    sizeY = maxY_ - minY_ + 1;
    sizeZ = maxZ_ - minZ_ + 1;
    pd.points = sizeX * sizeY * sizeZ;
    cd.cells  = (sizeX-1) * (sizeY-1) * (sizeZ-1);
}

void VTIWriter::setOrigin(double x_,double y_, double z_)
{
    originX = x_;
    originY = y_;
    originZ = z_;
}

void VTIWriter::setSpacing(double sx_,double sy_, double sz_)
{
    spaceX = sx_;
    spaceY = sy_;
    spaceZ = sz_;
}

void PVTIWriter::write()
{
    std::ofstream file;
    file.open(this -> filename);
    file << header();
    
    for (int i = 0; i< cellDataName.size(); i++)
    {
        if (i==0) file << "<PCellData>" << endl; 
        file << format_string( "<PDataArray Name=\"%s\" NumberOfComponents=\"%d\" type=\"%s\" />" ,cellDataName[i].c_str() ,cellDataComponents[i] ,cellDataType[i].c_str() ) << endl;
        if (i== cellDataName.size()-1) file << "</PCellData>" << endl; 
    }
    
    for (int i = 0; i< pointDataName.size(); i++)
    {
        if (i==0) file << "<PPointData>" << endl; 
        file << format_string( "<PDataArray Name=\"%s\" NumberOfComponents=\"%d\" type=\"%s\" />" ,pointDataName[i].c_str() ,pointDataComponents[i] ,pointDataType[i].c_str() ) << endl;
        if (i== pointDataName.size()-1) file << "</PPointData>" << endl; 
    }


    for (int i = 0; i< pieceFilename.size(); i++)
    {
        file << format_string("<Piece Extent= \"%d %d %d %d %d %d\" Source=\"%s\" />", pieceMinX[i], pieceMaxX[i], pieceMinY[i] , pieceMaxY[i] , pieceMinZ[i] , pieceMaxZ[i] , pieceFilename[i].c_str() ) << endl;
    }

    file << footer();
    file.close();
}

void VTIWriter::write()
{
    file.open(this -> filename);
    file << header();

    if (pd.sizeChild() > 0)
    {
        file << pd.header();
        for (int n = 0; n < pd.sizeChild() ; n++)
        {
            DataArray* array = (DataArray*) pd.getChild(n);
            file << *array;
        }
        file << pd.footer();
    }

    if (cd.sizeChild() > 0)
    {
        file << cd.header();
        for (int n = 0; n < cd.sizeChild() ; n++)
        {
            DataArray* array = (DataArray*) cd.getChild(n);
            file << *array;
        }
        file << cd.footer();
    }

    file << "</Piece>" << endl;
    file << "</ImageData>" << endl;

    file << ad;
    file << footer();
    file.close();
}

void PVTIWriter::addVTIWriter(VTIWriter& write)
{    
    pieceFilename.push_back( write.filename );
    
    pieceMaxX.push_back( write.pieceMaxX );
    pieceMaxY.push_back( write.pieceMaxY );
    pieceMaxZ.push_back( write.pieceMaxZ );
    pieceMinX.push_back( write.pieceMinX );
    pieceMinY.push_back( write.pieceMinY );
    pieceMinZ.push_back( write.pieceMinZ );
    
    if (pieceFilename.size() == 1)
    {
        for (int n = 0; n < write.pd.sizeChild() ; n++)
        {
            DataArray* data = (DataArray*) write.pd.getChild(n);
            pointDataName.push_back( data -> name );
            pointDataType.push_back( data -> type);
            pointDataComponents.push_back( data -> components );        
        }

        for (int n = 0; n < write.cd.sizeChild() ; n++)
        {
            DataArray* data = (DataArray*) write.cd.getChild(n);
            cellDataName.push_back( data -> name );
            cellDataType.push_back( data -> type);
            cellDataComponents.push_back( data -> components );        
        }

        originX = write.originX;
        originY = write.originX;
        originZ = write.originX;

        spaceX = write.spaceX;
        spaceY = write.spaceY;
        spaceZ = write.spaceZ;
               
        wholeMaxX = write.wholeMaxX;
        wholeMaxY = write.wholeMaxY;
        wholeMaxZ = write.wholeMaxZ;
        wholeMinX = write.wholeMinX;
        wholeMinY = write.wholeMinY;
        wholeMinZ = write.wholeMinZ;
    }
    else
    {
        wholeMaxX = max(wholeMaxX,write.wholeMaxX);
        wholeMaxY = max(wholeMaxY,write.wholeMaxY);
        wholeMaxZ = max(wholeMaxZ,write.wholeMaxZ);
        wholeMinX = min(wholeMinX,write.wholeMinX);
        wholeMinY = min(wholeMinY,write.wholeMinY);
        wholeMinZ = min(wholeMinZ,write.wholeMinZ);
    }
    

}

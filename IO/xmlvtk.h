/*
 * Copyright (c) 2025 Diogo Nardelli Siebert <diogo.siebert@ufsc.br>
 *
 * Licensed under either of
 *   - Apache License, Version 2.0 (https://www.apache.org/licenses/LICENSE-2.0)
 *   - GNU General Public License, Version 3.0 or later (https://www.gnu.org/licenses/gpl-3.0.html)
 *
 * SPDX-License-Identifier: (Apache-2.0 OR GPL-3.0-or-later)
 */

#ifndef __XMLVTK_H_INCLUDED__
#define __XMLVTK_H_INCLUDED__

#include <string>
#include <vector>
#include <fstream>
#include <cstdint>
#include <memory>

#define headerType u_int64_t


/**
 * @brief Encodes binary data into Base64 format.
 * @param input Pointer to the input data.
 * @param len Length of the input data.
 * @return Base64-encoded string.
 */
std::string spc_base64_encode(const unsigned char* input, size_t len);

/**
 * @class Element
 * @brief Represents a hierarchical XML-like element supporting compression and VTK serialization.
 */
class Element
{
    public:
        static bool compress; 
/** @brief Returns the parent element. */
        Element* getParent() { return parent ;}
/** @brief Returns the nth child element. *//** ... */
        Element* getChild(int n) { return this->child[n].get(); }
/** @brief Adds a new child element and sets its parent. */
    	void addChild(std::unique_ptr<Element> e) {
        	e -> parent = this;
	        child.push_back(std::move( e ) );
        }

/** @brief Returns the number of child elements. */
        int sizeChild( ) {  return child.size(); }
/** @brief Generates a header string for VTK output. */
        virtual std::string header() = 0;
/** @brief Generates a footer string for VTK output. */
        virtual std::string footer() = 0;
    private:
        Element* parent;
        std::vector< std::unique_ptr<Element> > child;
};

class PointData: public Element
{
    public:
        std::string header();
        std::string footer();
        unsigned int points;
};

class CellData: public Element
{
    public:
        std::string header();
        std::string footer();
        unsigned int cells;
};

class DataArray : public Element
{
    public:
        DataArray(const std::string&  name,const std::string&  type, const std::string&  format, int components, unsigned char* pointer, uint64_t size);
        unsigned char* pointer;
        std::string header();
        std::string footer();
        void write( std::ofstream& file );
        friend std::ostream& operator<<(std::ostream& os, DataArray& obj);
        uint64_t typeSize;
        uint64_t dataSize ;
        uint64_t offset;
        uint64_t points;
        std::string name;
        std::string format;
        std::string type;
        int components;
    private:
        int mode = 0;


};

class AppendData : public Element
{
    public:
        AppendData();
        std::string header();
        std::string footer();
        unsigned int addData(unsigned char* pointer, unsigned int size);
        friend std::ostream& operator<<(std::ostream& os, AppendData& obj);
    private:
        std::vector<unsigned char*> pointerList;
        std::vector<headerType> sizeList;
        headerType totalSize;
};

class VTIWriter : public Element
{
    public:
        VTIWriter(const std::string& filename);
        std::string footer();
        std::string header();

        void write();
        void setWholeExtent(int64_t  minX_,int64_t  minY_, int64_t  minZ_,int64_t  maxX_,int64_t  maxY_, int64_t  maxZ_);
        void setPiece(int64_t minX_,int64_t  minY_, int64_t  minZ_,int64_t  maxX_,int64_t  maxY_, int64_t  maxZ_);
        void setOrigin(double x_,double y_, double z_);
        void setSpacing(double sx_,double sy_, double sz_);
        void addPointData(const std::string&  name,const std::string&  type, const std::string&  format, int components, unsigned char* pointer);
        void addCellData(const std::string&  name,const std::string&  type, const std::string&  format, int components, unsigned char* pointer);
        void setCompress();

        std::string filename;
        std::string vtkVersion;               /*!< The Vtk File Format Version of the file */
        std::string fileTitle;                /*!< The title of the file (do not confuse with the name of the file) */
        std::string dataSetType;              /*!< The type of geometry (grid) that data is associeted to (STRUCTURED GRID for LBM applications) */
        std::string headerTypeName;   

        int64_t   wholeMinX,  wholeMinY,  wholeMinZ;
        int64_t   wholeMaxX,  wholeMaxY,  wholeMaxZ;
        int64_t   pieceMinX,  pieceMinY,  pieceMinZ;
        int64_t   pieceMaxX,  pieceMaxY,  pieceMaxZ;
        int64_t   sizeX,  sizeY,  sizeZ;         /*!< Lenght in pixels of the image in each axis */

        /*!< Lenght in pixels of the image in each axis */
        double spaceX, spaceY, spaceZ;          /*!< Ratio of the different axis */
        double originX,originY,originZ;          /*!< Position of the origin of the image */

        std::string byteOrder;       /*!< Position of the origin of the image */
        std::ofstream file;
	    CellData cd;
        PointData pd;
        AppendData ad;
};

class PVTIWriter : public Element
{
        public:
            PVTIWriter(const std::string& filename);
            std::string filename;
            std::string footer();
            std::string header();

            double spaceX, spaceY, spaceZ;      
            double originX,originY,originZ;     
            
            int pieceCounter;

            int64_t   wholeMinX,  wholeMinY,  wholeMinZ;
            int64_t   wholeMaxX,  wholeMaxY,  wholeMaxZ;

            std::vector<int64_t> pieceMinX,  pieceMinY,  pieceMinZ;
            std::vector<int64_t> pieceMaxX,  pieceMaxY,  pieceMaxZ;
            std::vector<std::string> pieceFilename;

            std::vector<int> cellDataComponents;
            std::vector<std::string> cellDataName;
            std::vector<std::string> cellDataType;

            std::vector<int> pointDataComponents;
            std::vector<std::string> pointDataName;
            std::vector<std::string> pointDataType;

            void addVTIWriter(VTIWriter& write);
            void write();
};


#endif

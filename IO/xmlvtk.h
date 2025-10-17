#ifndef __XMLVTK_H_INCLUDED__
#define __XMLVTK_H_INCLUDED__

#define BASE_64 0
#define APPEND_RAW_DATA 1

#include <string>
#include <vector>
#include <fstream>
#include <cstdint>

#define headerType u_int64_t

static char b64table[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

unsigned char *spc_base64_encode(unsigned char *input, size_t len);

class Element
{
    public:
        static bool compress; 
        Element* getParent() { return parent ;}
        Element* getChild(int n) { return this->child[n] ; }
        void addChild( Element* e) { child.push_back(e); e->parent = this; }
        int sizeChild( ) {  return child.size(); }
        virtual std::string header() = 0;
        virtual std::string footer() = 0;
    private:
        Element* parent;
        std::vector<Element*> child;
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
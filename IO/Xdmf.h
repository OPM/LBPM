#include "IO/HDF5_IO.h"
#include "common/Array.h"
#include "common/MPI.h"
#include "common/Utilities.h"

#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <vector>


// Helper class to write/read XDMF files
class Xdmf
{
public:
    enum class TopologyType {
        Null = 0,
        Polyvertex,
        Polyline,
        Polygon,
        Triangle,
        Quadrilateral,
        Tetrahedron,
        Pyramid,
        Wedge,
        Hexahedron,
        Edge_3,
        Triangle_6,
        Quadrilateral_8,
        Tetrahedron_10,
        Pyramid_13,
        Wedge_15,
        Hexahedron_20,
        Mixed,
        CurvilinearMesh2D,
        CurvilinearMesh3D,
        RectangularMesh2D,
        RectangularMesh3D,
        UniformMesh2D,
        UniformMesh3D,
    };
    enum class DataType { Null = 0, Char, Int32, Int64, Uint32, Uint64, Float, Double };
    enum class RankType { Null = 0, Scalar, Vector, Tensor, Tensor6, Matrix, GlobalID };
    enum class Center { Null = 0, Node, Edge, Face, Cell, Grid, Other };

    struct VarData {
        std::string name;  // Variable name
        ArraySize size;    // Size of variable
        RankType rankType; // Rank order of data
        Center center;     // Variable centering
        std::string data;  // Variable data
    };

    struct MeshData {
        std::string name;          // Name of mesh
        TopologyType type;         // Type of mesh
        ArraySize size;            // Size of mesh (meaning depends on mesh type)
        double range[6];           // Range of the mesh (only used for UniformMesh2D/UniformMesh3D)
        std::string x;             // x coordinates (or xy/xyz coordinates)
        std::string y;             // y coordinates
        std::string z;             // z coordinates
        std::string dofMap;        // mesh connectivity
        std::vector<VarData> vars; // Variables
        MeshData() : type( TopologyType::Null ), range{ 0 } {}
        //! Add a variable
        void addVariable( const std::string &meshName, const std::string &varName,
            ArraySize varSize, RankType rank, Center center, const std::string &varData );
    };


public:
    //! Add a Point mesh
    static MeshData createPointMesh( const std::string &name, uint8_t NDIM, size_t N,
        const std::string &x, const std::string &y = "", const std::string &z = "" );

    /*!
     * @brief  Add a uniform mesh
     * @details  This function adds a uniform rectangular mesh
     * @param[in] name          The name of the mesh
     * @param[in] range         The range of the mesh [ x_min, x_max, y_min, y_max, z_min, z_max ]
     * @param[in] size          The number of cells in the mesh
     */
    static MeshData createUniformMesh(
        const std::string &name, const std::vector<double> &range, ArraySize size );

    /*!
     * @brief  Add a Curvilinear mesh
     * @details  This function adds a curvilinear mesh
     * @param[in] name          The name of the mesh
     * @param[in] size          The number of cells in the mesh
     * @param[in] x             The x coordinates or the xy/xyz coordinates
     * @param[in] y             The y coordinates (may be null)
     * @param[in] z             The z coordinates (may be null)
     */
    static MeshData createCurvilinearMesh( const std::string &name, ArraySize size,
        const std::string &x, const std::string &y, const std::string &z = "" );

    /*!
     * @brief  Add an unstructured mesh
     * @details  This function adds an unstructerd mesh to the class to write.
     *    The mesh may be one of several unsupported unstructured mesh types.
     *    This function does not support mixed elements.
     * @param[in] name          The name of the mesh
     * @param[in] NDIM          The number of physical dimensions
     * @param[in] type          The element type
     * @param[in] NumElements   The number of elements
     * @param[in] dofMap        The connectivity information (type x NumElements)
     * @param[in] NumNodes      The number of nodes
     * @param[in] x             The x coordinates or the xy/xyz coordinates
     * @param[in] y             The y coordinates (may be null)
     * @param[in] z             The z coordinates (may be null)
     */
    static MeshData createUnstructuredMesh( const std::string &name, uint8_t NDIM,
        TopologyType type, size_t NumElements, const std::string &dofMap, size_t NumNodes,
        const std::string &x, const std::string &y = "", const std::string &z = "" );


public:
    //! Add a sub-domain
    void addMesh( const std::string &meshName, const MeshData &domain );

    //! Gather all data to rank 0
    void gather( const Utilities::MPI &comm );

    //! Write the xml file
    void write( const std::string &filename ) const;


private:
    std::map<std::string, std::vector<MeshData>> d_meshData;
};

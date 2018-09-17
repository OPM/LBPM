#include <vector>
#include "analysis/pmmc.h"

/* 
Doubly-connected edge list (DECL) 
*/

// Vertex structure
class Vertex{
public:
	Vertex() { d_data.resize(12); }
	~Vertex() = default;
    Vertex( const Vertex& ) = delete;
    Vertex operator=( const Vertex& ) = delete;

    // Add/assign a point
	inline void add( const Point& P ) { d_data.push_back( P ); }
	inline void assign( int idx, const Point& P ) { d_data[idx] = P; }

    // Get a point
	inline Point& coords( int idx ) { return d_data[idx]; }
	inline const Point& coords( int idx ) const { return d_data[idx]; }

	int IncidentEdge();

    // Return the number of points
	inline int size() const { return d_data.size(); }

private:
	std::vector<Point> d_data;
};


// Halfedge structure
// Face
class Halfedge{
public:
	Halfedge() = default;
	~Halfedge() = default;
    Halfedge( const Halfedge& ) = delete;
    Halfedge operator=( const Halfedge& ) = delete;

	inline int v1(int edge) const { return d_data[edge][0]; }
	inline int v2(int edge) const { return d_data[edge][1]; }
	inline int face(int edge) const { return d_data[edge][2]; }
	inline int twin(int edge) const { return d_data[edge][3]; }
	inline int prev(int edge) const { return d_data[edge][4]; }
	inline int next(int edge) const { return d_data[edge][5]; }

	inline int size() const { return d_data.size(); }
    inline void resize( int N ) { d_data.resize( N ); }

    inline int& data( int i, int j ) { return d_data[j][i]; }
    inline const int& data( int i, int j ) const { return d_data[j][i]; }

private:
	std::vector<std::array<int,6>> d_data;
};

// DECL
class DECL{
public:
	DECL();
	~DECL();
	
	int face();
	Vertex vertex;
	Halfedge halfedge;
	void LocalIsosurface(const DoubleArray& A, double value, int i, int j, int k);
	int Face(int index);
	
	double origin(int edge);
	double EdgeAngle(int edge);
	Point TriNormal(int edge);
	int TriangleCount;
	int VertexCount;

private:
	std::vector<int> FaceData;
};

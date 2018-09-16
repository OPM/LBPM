#include <vector>
#include "analysis/pmmc.h"

/* 
Doubly-connected edge list (DECL) 
*/

// Vertex structure
class Vertex{
public:
	Vertex();
	~Vertex();
	void add(Point P);
	void assign( int idx, Point P);
	int size();
	Point coords(int idx);
	int IncidentEdge();
private:
	std::vector<double> vertex_data;
	int size_;
};

// Halfedge structure
// Face
class Halfedge{
public:
	Halfedge();
	~Halfedge();

	int v1(int edge);
	int v2(int edge);
	int twin(int edge);
	int face(int edge);
	int next(int edge);
	int prev(int edge);
	int size();

	Array<int> data;
private:
	int size_;
};

// DECL
class DECL{
public:
	DECL();
	~DECL();
	
	int face();
	Vertex vertex;
	Halfedge halfedge;
	void LocalIsosurface(const DoubleArray A, double value, int i, int j, int k);
	int Face(int index);
	
	double origin(int edge);
	double EdgeAngle(int edge);
	Point TriNormal(int edge);
	int TriangleCount;
	int VertexCount;	

private:
	Array <int> FaceData;

	Point VertexList[12];
	Point NewVertexList[12];
	int LocalRemap[12];
	double CubeValues[8];	

	DTMutableList<Point> cellvertices;// = DTMutableList<Point>(20);
	IntArray Triangles;// = IntArray(3,20);

};

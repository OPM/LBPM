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
	void assign(unsigned long int idx, Point P);
	unsigned long int size();
	Point coords(unsigned long int idx);
	unsigned long int IncidentEdge();
private:
	std::vector<double> vertex_data;
	unsigned long int size_;
};

// Halfedge structure
// Face
class Halfedge{
public:
	Halfedge();
	~Halfedge();

	unsigned long int v1(unsigned long int edge);
	unsigned long int v2(unsigned long int edge);
	unsigned long int twin(unsigned long int edge);
	unsigned long int face(unsigned long int edge);
	unsigned long int next(unsigned long int edge);
	unsigned long int prev(unsigned long int edge);
	unsigned long int size();

	Array<unsigned long int> data;
private:
	unsigned long int size_;
};

// DECL
class DECL{
public:
	DECL();
	~DECL();
	
	unsigned long int face();
	Vertex vertex;
	Halfedge halfedge;
	void LocalIsosurface(const DoubleArray A, double value, int i, int j, int k);
	
	double origin(int edge);
	double EdgeAngle(int edge);
	Point TriNormal(int edge);
	unsigned long int TriangleCount;
	unsigned long int VertexCount;	

private:
	unsigned long int *face_data;


};

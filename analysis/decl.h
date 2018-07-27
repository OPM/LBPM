#include <vector>

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
	Point coords(unsigned long int idx);
	unsigned long int IncidentEdge();
	unsigned long int count;
private:
	std::vector<double> vertex_data;
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
	
private:
	Array<unsigned long int> HalfEdge;
};

// DECL
class DECL{
public:
	DECL();
	~DECL();
	
	unsigned long int face();
	Vertex vertex;
	Halfedge halfedge;
	void AddCube(); // need a function to add new faces based on marching cubes surface
	
	double origin(int edge);
	double EdgeAngle(int edge);
	Point TriNormal(int edge);
	
private:
	unsigned long int *face_data;
};

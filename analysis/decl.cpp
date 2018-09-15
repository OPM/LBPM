#include "analysis/decl.h"

/* 
Double connected edge list (DECL) 
*/

Vertex::Vertex(){
	size_ = 0;
	vertex_data.resize(36);
}

Vertex::~Vertex(){

}

void Vertex::add(Point P){
	vertex_data.push_back(P.x);
	vertex_data.push_back(P.y);
	vertex_data.push_back(P.z);
	size_++;
}

void Vertex::assign(int idx, Point P){
	vertex_data[3*idx] = P.x;
	vertex_data[3*idx+1] = P.y;
	vertex_data[3*idx+2] = P.z;
}

int Vertex::size(){
	return size_;
}

Point Vertex::coords(int idx){
	Point P;
	P.x = vertex_data[3*idx];
	P.y = vertex_data[3*idx+1];
	P.z = vertex_data[3*idx+2];
	return P;
}


Halfedge::Halfedge(){
	size_=0;
}

Halfedge::~Halfedge(){
	
}
int Halfedge::v1(int edge){
	return data(0,edge);
}

int Halfedge::v2(int edge){
	return data(1,edge);
}

int Halfedge::face(int edge){
	return data(2,edge);
}

int Halfedge::twin(int edge){
	return data(3,edge);
}

int Halfedge::prev(int edge){
	return data(4,edge);
}

int Halfedge::next(int edge){
	return data(5,edge);
}

int Halfedge::size(){
	return size_;
}


DECL::DECL(){
}

DECL::~DECL(){
	TriangleCount=0;
	VertexCount=0;
	
}

int DECL::Face(int index){
	return FaceData(index);
}

void DECL::LocalIsosurface(const DoubleArray A, double value, const int i, const int j, const int k){
	Point P,Q;
	Point PlaceHolder;
	Point C0,C1,C2,C3,C4,C5,C6,C7;

	int CubeIndex;
	int nTris = 0;
	int nVert =0;

	Point VertexList[12];
	Point NewVertexList[12];
	int LocalRemap[12];

	DTMutableList<Point> cellvertices = DTMutableList<Point>(20);
	IntArray Triangles = IntArray(3,20);

	// Values from array 'A' at the cube corners
	double CubeValues[8];

	// Points corresponding to cube corners
	C0.x = 0.0; C0.y = 0.0; C0.z = 0.0;
	C1.x = 1.0; C1.y = 0.0; C1.z = 0.0;
	C2.x = 1.0; C2.y = 1.0; C2.z = 0.0;
	C3.x = 0.0; C3.y = 1.0; C3.z = 0.0;
	C4.x = 0.0; C4.y = 0.0; C4.z = 1.0;
	C5.x = 1.0; C5.y = 0.0; C5.z = 1.0;
	C6.x = 1.0; C6.y = 1.0; C6.z = 1.0;
	C7.x = 0.0; C7.y = 1.0; C7.z = 1.0;							

	CubeValues[0] = A(i,j,k) - value;
	CubeValues[1] = A(i+1,j,k) - value;
	CubeValues[2] = A(i+1,j+1,k) - value;
	CubeValues[3] = A(i,j+1,k) - value;
	CubeValues[4] = A(i,j,k+1) - value;
	CubeValues[5] = A(i+1,j,k+1) - value;
	CubeValues[6] = A(i+1,j+1,k+1) - value;
	CubeValues[7] = A(i,j+1,k+1) -value;
	//printf("Set cube values: %i, %i, %i \n",i,j,k);
	
	//Determine the index into the edge table which
	//tells us which vertices are inside of the surface
	CubeIndex = 0;
	if (CubeValues[0] < 0.0f) CubeIndex |= 1;
	if (CubeValues[1] < 0.0f) CubeIndex |= 2;
	if (CubeValues[2] < 0.0f) CubeIndex |= 4;
	if (CubeValues[3] < 0.0f) CubeIndex |= 8;
	if (CubeValues[4] < 0.0f) CubeIndex |= 16;
	if (CubeValues[5] < 0.0f) CubeIndex |= 32;
	if (CubeValues[6] < 0.0f) CubeIndex |= 64;
	if (CubeValues[7] < 0.0f) CubeIndex |= 128;

	//Find the vertices where the surface intersects the cube
	if (edgeTable[CubeIndex] & 1){
		P = VertexInterp(C0,C1,CubeValues[0],CubeValues[1]);
		VertexList[0] = P;
		Q = C0;
	}
	if (edgeTable[CubeIndex] & 2){
		P = VertexInterp(C1,C2,CubeValues[1],CubeValues[2]);
		VertexList[1] = P;
		Q = C1;
	}
	if (edgeTable[CubeIndex] & 4){
		P = VertexInterp(C2,C3,CubeValues[2],CubeValues[3]);
		VertexList[2] =	P;
		Q = C2;
	}
	if (edgeTable[CubeIndex] & 8){
		P = VertexInterp(C3,C0,CubeValues[3],CubeValues[0]);
		VertexList[3] =	P;
		Q = C3;
	}
	if (edgeTable[CubeIndex] & 16){
		P = VertexInterp(C4,C5,CubeValues[4],CubeValues[5]);
		VertexList[4] =	P;
		Q = C4;
	}
	if (edgeTable[CubeIndex] & 32){
		P = VertexInterp(C5,C6,CubeValues[5],CubeValues[6]);
		VertexList[5] = P;
		Q = C5;
	}
	if (edgeTable[CubeIndex] & 64){ 
		P = VertexInterp(C6,C7,CubeValues[6],CubeValues[7]);
		VertexList[6] = P;
		Q = C6;
	}
	if (edgeTable[CubeIndex] & 128){
		P = VertexInterp(C7,C4,CubeValues[7],CubeValues[4]);
		VertexList[7] =	P;
		Q = C7;
	}
	if (edgeTable[CubeIndex] & 256){
		P = VertexInterp(C0,C4,CubeValues[0],CubeValues[4]);
		VertexList[8] =	P;
		Q = C0;
	}
	if (edgeTable[CubeIndex] & 512){
		P = VertexInterp(C1,C5,CubeValues[1],CubeValues[5]);
		VertexList[9] =	P;
		Q = C1;
	}
	if (edgeTable[CubeIndex] & 1024){
		P = VertexInterp(C2,C6,CubeValues[2],CubeValues[6]);
		VertexList[10] = P;
		Q = C2;
	}
	if (edgeTable[CubeIndex] & 2048){
		P = VertexInterp(C3,C7,CubeValues[3],CubeValues[7]);
		VertexList[11] = P;
		Q = C3;
	}

	VertexCount=0;
	for (int idx=0;idx<12;idx++)
		LocalRemap[idx] = -1;

	for (int idx=0;triTable[CubeIndex][idx]!=-1;idx++)
	{
		if(LocalRemap[triTable[CubeIndex][idx]] == -1)
		{
			NewVertexList[VertexCount] = VertexList[triTable[CubeIndex][idx]];
			LocalRemap[triTable[CubeIndex][idx]] = VertexCount;
			VertexCount++;
		}
	}
	
	//printf("Found %i vertices \n",VertexCount);

	for (int idx=0;idx<VertexCount;idx++) {
		P = NewVertexList[idx];
		//P.x  += i;
		//P.y  += j;
		//P.z  += k;
		cellvertices(idx) = P;
	}
	nVert = VertexCount;

	TriangleCount = 0;
	for (int idx=0;triTable[CubeIndex][idx]!=-1;idx+=3) {
		Triangles(0,TriangleCount) = LocalRemap[triTable[CubeIndex][idx+0]];
		Triangles(1,TriangleCount) = LocalRemap[triTable[CubeIndex][idx+1]];
		Triangles(2,TriangleCount) = LocalRemap[triTable[CubeIndex][idx+2]];
		TriangleCount++;
	}
	nTris = TriangleCount;

	// Now add the local values to the DECL data structure
	if (nTris>0){
		FaceData.resize(TriangleCount);
	    //printf("Construct halfedge structure... \n");
	    //printf("   Construct %i triangles \n",nTris);
		halfedge.data.resize(6,nTris*3);
		int idx_edge=0;
		for (int idx=0; idx<TriangleCount; idx++){
			int V1 = Triangles(0,idx);
			int V2 = Triangles(1,idx);
			int V3 = Triangles(2,idx);
			FaceData(idx) = idx_edge;
			// first edge: V1->V2 
			halfedge.data(0,idx_edge) = V1;  // first vertex
			halfedge.data(1,idx_edge) = V2;  // second vertex
			halfedge.data(2,idx_edge) = idx; // triangle
			halfedge.data(3,idx_edge) = -1;  // twin
			halfedge.data(4,idx_edge) = idx_edge+2;  // previous edge
			halfedge.data(5,idx_edge) = idx_edge+1;  // next edge
			idx_edge++;
			// second edge: V2->V3
			halfedge.data(0,idx_edge) = V2;  // first vertex
			halfedge.data(1,idx_edge) = V3;  // second vertex
			halfedge.data(2,idx_edge) = idx; // triangle
			halfedge.data(3,idx_edge) = -1;  // twin
			halfedge.data(4,idx_edge) = idx_edge-1;  // previous edge
			halfedge.data(5,idx_edge) = idx_edge+1;  // next edge
			idx_edge++;
			// third edge: V3->V1
			halfedge.data(0,idx_edge) = V3;  // first vertex
			halfedge.data(1,idx_edge) = V1;  // second vertex
			halfedge.data(2,idx_edge) = idx; // triangle
			halfedge.data(3,idx_edge) = -1;  // twin
			halfedge.data(4,idx_edge) = idx_edge-1;  // previous edge
			halfedge.data(5,idx_edge) = idx_edge-2;  // next edge
			idx_edge++;
			//printf("   ***tri %i ***edge %i *** \n",idx, idx_edge);
		}
		//printf("  parsing halfedge structure\n");
		int EdgeCount=idx_edge;
		for (int idx=0; idx<EdgeCount; idx++){
			int V1=halfedge.data(0,idx);
			int V2=halfedge.data(1,idx);
			// Find all the twins within the cube
			for (int jdx=0; jdx<EdgeCount; jdx++){
				if (halfedge.data(1,jdx) == V1 && halfedge.data(0,jdx) == V2){
					// this is the pair
					halfedge.data(3,idx) = jdx;
					halfedge.data(3,jdx) = idx;
				}
				if (halfedge.data(1,jdx) == V2 && halfedge.data(0,jdx) == V1 && !(idx==jdx)){
					std::printf("WARNING: half edges with identical orientation! \n");
				}
			}	
			// Use "ghost" twins if edge is on a cube face
			P = cellvertices(V1);
			Q = cellvertices(V2);
			if (P.x == 0.0 && Q.x == 0.0) halfedge.data(3,idx) = -1;  // ghost twin for x=0 face
			if (P.x == 1.0 && Q.x == 1.0) halfedge.data(3,idx) = -4;  // ghost twin for x=1 face
			if (P.y == 0.0 && Q.y == 0.0) halfedge.data(3,idx) = -2;  // ghost twin for y=0 face
			if (P.y == 1.0 && Q.y == 1.0) halfedge.data(3,idx) = -5;  // ghost twin for y=1 face
			if (P.z == 0.0 && Q.z == 0.0) halfedge.data(3,idx) = -3;  // ghost twin for z=0 face
			if (P.z == 1.0 && Q.z == 1.0) halfedge.data(3,idx) = -6;  // ghost twin for z=1 face
		}
	}

	// Map vertices to global coordinates
	for (int idx=0;idx<VertexCount;idx++) {
		P = cellvertices(idx);
		P.x  += i;
		P.y  += j;
		P.z  += k;
		vertex.assign(idx,P);
	}
}

Point DECL::TriNormal(int edge)
{
	Point P,Q,R;
	Point U,V,W;
	double ux,uy,uz,vx,vy,vz;
	double nx,ny,nz,len;
	// at cube faces define outward normal to cube
	if (edge == -1){
		W.x = -1.0; W.y = 0.0; W.z = 0.0; // x cube face
	}
	else if (edge == -2){
		W.x = 0.0; W.y = -1.0; W.z = 0.0; // y cube face
	}
	else if (edge == -3){
		W.x = 0.0; W.y = 0.0; W.z = -1.0; // z cube face
	}
	else if (edge == -4){
		W.x = 1.0; W.y = 0.0; W.z = 0.0; // x cube face
	}
	else if (edge == -5){
		W.x = 0.0; W.y = 1.0; W.z = 0.0; // y cube face
	}
	else if (edge == -6){
		W.x = 0.0; W.y = 0.0; W.z = 1.0; // z cube face
	}
	else{
		// vertices for triange
		int e2 = halfedge.next(edge);
		int e3 = halfedge.next(e2);
		P=vertex.coords(halfedge.v1(edge));
		Q=vertex.coords(halfedge.v1(e2));
		R=vertex.coords(halfedge.v1(e3));
		// edge vectors
		U = Q-P;
		V = R-Q;
		// normal vector
		nx = U.y*V.z - U.z*V.y;
		ny = U.z*V.x - U.x*V.z;
		nz = U.x*V.y - U.y*V.x;
		len = sqrt(nx*nx+ny*ny+nz*nz);
		W.x = nx/len; W.y = ny/len; W.z = nz/len;
	}
	return W;
}

double DECL::EdgeAngle(int edge)
{
	double angle;
	double dotprod,length,hypotenuse;
	Point P,Q,R; // triangle vertices
	Point U,V,W; // normal vectors
	int e2 = halfedge.next(edge);
	int e3 = halfedge.next(e2);
	P=vertex.coords(halfedge.v1(edge));
	Q=vertex.coords(halfedge.v1(e2));
	R=vertex.coords(halfedge.v1(e3));
	U = TriNormal(edge);
	V = TriNormal(halfedge.twin(edge));
	if (halfedge.twin(edge) < 0 ){
		// compute edge normal in plane of cube face
		W = P - Q; // edge tangent vector
		length = sqrt(W.x*W.x+W.y*W.y+W.z*W.z);
		W.x /= length;
		W.y /= length;
		W.z /= length;
		// edge normal within the plane of the cube face
		double nx = W.y*V.z - W.z*V.y;
		double ny = W.z*V.x - W.x*V.z;
		double nz = W.x*V.y - W.y*V.x;
		double len = sqrt(nx*nx+ny*ny+nz*nz);
		// new value for V is this normal vector
		V.x = nx/len; V.y = ny/len; V.z = nz/len;
		dotprod = U.x*V.x + U.y*V.y + U.z*V.z;
		if (dotprod > 1.f) dotprod=1.f;
		if (dotprod < -1.f) dotprod=-1.f;
		angle = acos(dotprod);
		/* project onto plane of cube face also works
		W = U - dotprod*V;
		length = sqrt(W.x*W.x+W.y*W.y+W.z*W.z); // for normalization
		dotprod = (U.x*W.x + U.y*W.y + U.z*W.z)/length;
		if (dotprod > 1.f) dotprod=1.f;
		if (dotprod < -1.f) dotprod=-1.f;
		angle = acos(dotprod);
		*/
	}
	else{
		dotprod=U.x*V.x + U.y*V.y + U.z*V.z;
		if (dotprod > 1.f) dotprod=1.f;
		if (dotprod < -1.f) dotprod=-1.f;
		angle = 0.5*acos(dotprod);
	}
	// determine if angle is concave or convex based on edge normal
	W = 0.5*(P+Q)-R; // vector that lies in plane of triangle
	//hypotenuse = sqrt(W.x*W.x+W.y*W.y+W.z*W.z + V.x*V.x+V.y*V.y+V.z*V.z); // hypotenuse of right triangle
	//length = sqrt((W.x+V.x)*(W.x+V.x) + (W.y+V.y)*(W.y+V.y) + (W.z+V.z)*(W.z+V.z));
	//if (length > hypotenuse){
	if (W.x*V.x + W.y*V.y + W.z*V.z < 0.f){
		// concave
		angle = -angle;
	}
	//printf("     %f, %f (Edge=%i, twin=%i)\n       U={%f, %f, %f}, V={%f, %f, %f}\n",angle,dotprod,edge,halfedge.twin(edge),U.x,U.y,U.z,V.x,V.y,V.z);
	return angle;
}

void Isosurface(DoubleArray &A, const double &v)
{
	Point P,Q;
	Point PlaceHolder;
	double temp;
	Point C0,C1,C2,C3,C4,C5,C6,C7;

	int TriangleCount;
	int VertexCount;
	int CubeIndex;
	int nTris, nVert;

	Point VertexList[12];
	Point NewVertexList[12];
	int LocalRemap[12];

	DTMutableList<Point> cellvertices = DTMutableList<Point>(20);
	IntArray Triangles = IntArray(3,20);

	// Values from array 'A' at the cube corners
	double CubeValues[8];

	int Nx = A.size(0);
	int Ny = A.size(1);
	int Nz = A.size(2);

	// Points corresponding to cube corners
	C0.x = 0.0; C0.y = 0.0; C0.z = 0.0;
	C1.x = 1.0; C1.y = 0.0; C1.z = 0.0;
	C2.x = 1.0; C2.y = 1.0; C2.z = 0.0;
	C3.x = 0.0; C3.y = 1.0; C3.z = 0.0;
	C4.x = 0.0; C4.y = 0.0; C4.z = 1.0;
	C5.x = 1.0; C5.y = 0.0; C5.z = 1.0;
	C6.x = 1.0; C6.y = 1.0; C6.z = 1.0;
	C7.x = 0.0; C7.y = 1.0; C7.z = 1.0;							

	for (int k=1; k<Nz-1; k++){		
		for (int j=1; j<Ny-1; j++){
			for (int i=1; i<Nx-1; i++){
				// Set the corner values for this cube
				CubeValues[0] = A(i,j,k);
				CubeValues[1] = A(i+1,j,k);
				CubeValues[2] = A(i+1,j+1,k);
				CubeValues[3] = A(i,j+1,k);
				CubeValues[4] = A(i,j,k+1);
				CubeValues[5] = A(i+1,j,k+1);
				CubeValues[6] = A(i+1,j+1,k+1);
				CubeValues[7] = A(i,j+1,k+1);

				//Determine the index into the edge table which
				//tells us which vertices are inside of the surface
				CubeIndex = 0;
				if (CubeValues[0] < 0.0f) CubeIndex |= 1;
				if (CubeValues[1] < 0.0f) CubeIndex |= 2;
				if (CubeValues[2] < 0.0f) CubeIndex |= 4;
				if (CubeValues[3] < 0.0f) CubeIndex |= 8;
				if (CubeValues[4] < 0.0f) CubeIndex |= 16;
				if (CubeValues[5] < 0.0f) CubeIndex |= 32;
				if (CubeValues[6] < 0.0f) CubeIndex |= 64;
				if (CubeValues[7] < 0.0f) CubeIndex |= 128;

				//Find the vertices where the surface intersects the cube
				if (edgeTable[CubeIndex] & 1){
					P = VertexInterp(C0,C1,CubeValues[0],CubeValues[1]);
					VertexList[0] = P;
					Q = C0;
				}
				if (edgeTable[CubeIndex] & 2){
					P = VertexInterp(C1,C2,CubeValues[1],CubeValues[2]);
					VertexList[1] = P;
					Q = C1;
				}
				if (edgeTable[CubeIndex] & 4){
					P = VertexInterp(C2,C3,CubeValues[2],CubeValues[3]);
					VertexList[2] =	P;
					Q = C2;
				}
				if (edgeTable[CubeIndex] & 8){
					P = VertexInterp(C3,C0,CubeValues[3],CubeValues[0]);
					VertexList[3] =	P;
					Q = C3;
				}
				if (edgeTable[CubeIndex] & 16){
					P = VertexInterp(C4,C5,CubeValues[4],CubeValues[5]);
					VertexList[4] =	P;
					Q = C4;
				}
				if (edgeTable[CubeIndex] & 32){
					P = VertexInterp(C5,C6,CubeValues[5],CubeValues[6]);
					VertexList[5] = P;
					Q = C5;
				}
				if (edgeTable[CubeIndex] & 64){ 
					P = VertexInterp(C6,C7,CubeValues[6],CubeValues[7]);
					VertexList[6] = P;
					Q = C6;
				}
				if (edgeTable[CubeIndex] & 128){
					P = VertexInterp(C7,C4,CubeValues[7],CubeValues[4]);
					VertexList[7] =	P;
					Q = C7;
				}
				if (edgeTable[CubeIndex] & 256){
					P = VertexInterp(C0,C4,CubeValues[0],CubeValues[4]);
					VertexList[8] =	P;
					Q = C0;
				}
				if (edgeTable[CubeIndex] & 512){
					P = VertexInterp(C1,C5,CubeValues[1],CubeValues[5]);
					VertexList[9] =	P;
					Q = C1;
				}
				if (edgeTable[CubeIndex] & 1024){
					P = VertexInterp(C2,C6,CubeValues[2],CubeValues[6]);
					VertexList[10] = P;
					Q = C2;
				}
				if (edgeTable[CubeIndex] & 2048){
					P = VertexInterp(C3,C7,CubeValues[3],CubeValues[7]);
					VertexList[11] = P;
					Q = C3;
				}

				VertexCount=0;
				for (int idx=0;idx<12;idx++)
					LocalRemap[idx] = -1;

				for (int idx=0;triTable[CubeIndex][idx]!=-1;idx++)
				{
					if(LocalRemap[triTable[CubeIndex][idx]] == -1)
					{
						NewVertexList[VertexCount] = VertexList[triTable[CubeIndex][idx]];
						LocalRemap[triTable[CubeIndex][idx]] = VertexCount;
						VertexCount++;
					}
				}

				for (int idx=0;idx<VertexCount;idx++) {
					P = NewVertexList[idx];
					//P.x  += i;
					//P.y  += j;
					//P.z  += k;
					cellvertices(idx) = P;
				}
				nVert = VertexCount;

				TriangleCount = 0;
				for (int idx=0;triTable[CubeIndex][idx]!=-1;idx+=3) {
					Triangles(0,TriangleCount) = LocalRemap[triTable[CubeIndex][idx+0]];
					Triangles(1,TriangleCount) = LocalRemap[triTable[CubeIndex][idx+1]];
					Triangles(2,TriangleCount) = LocalRemap[triTable[CubeIndex][idx+2]];
					TriangleCount++;
				}
				nTris = TriangleCount;

				// Now add the local values to the DECL data structure
				IntArray HalfEdge(6,nTris*3);
				DoubleArray EdgeAngles(nTris*3);
				int idx_edge=0;
				for (int idx=0; idx<TriangleCount; idx++){
					int V1 = Triangles(0,idx);
					int V2 = Triangles(1,idx);
					int V3 = Triangles(2,idx);
					// first edge: V1->V2 
					HalfEdge(0,idx_edge) = V1;  // first vertex
					HalfEdge(1,idx_edge) = V2;  // second vertex
					HalfEdge(2,idx_edge) = idx; // triangle
					HalfEdge(3,idx_edge) = -1;  // twin
					HalfEdge(4,idx_edge) = idx_edge+2;  // previous edge
					HalfEdge(5,idx_edge) = idx_edge+1;  // next edge
					idx_edge++;
					// second edge: V2->V3
					HalfEdge(0,idx_edge) = V2;  // first vertex
					HalfEdge(1,idx_edge) = V3;  // second vertex
					HalfEdge(2,idx_edge) = idx; // triangle
					HalfEdge(3,idx_edge) = -1;  // twin
					HalfEdge(4,idx_edge) = idx_edge-1;  // previous edge
					HalfEdge(5,idx_edge) = idx_edge+1;  // next edge
					idx_edge++;
					// third edge: V3->V1
					HalfEdge(0,idx_edge) = V3;  // first vertex
					HalfEdge(1,idx_edge) = V1;  // second vertex
					HalfEdge(2,idx_edge) = idx; // triangle
					HalfEdge(3,idx_edge) = -1;  // twin
					HalfEdge(4,idx_edge) = idx_edge-1;  // previous edge
					HalfEdge(5,idx_edge) = idx_edge-2;  // next edge
					idx_edge++;
				}
				int EdgeCount=idx_edge;
				for (int idx=0; idx<EdgeCount; idx++){
					int V1=HalfEdge(0,idx);
					int V2=HalfEdge(1,idx);
					// Find all the twins within the cube
					for (int jdx=0; idx<EdgeCount; jdx++){
						if (HalfEdge(1,jdx) == V1 && HalfEdge(0,jdx) == V2){
							// this is the pair
							HalfEdge(3,idx) = jdx;
							HalfEdge(3,jdx) = idx;
						}
						if (HalfEdge(1,jdx) == V2 && HalfEdge(0,jdx) == V1 && !(idx==jdx)){
							std::printf("WARNING: half edges with identical orientation! \n");
						}
					}	
					// Use "ghost" twins if edge is on a cube face
					P = cellvertices(V1);
					Q = cellvertices(V2);
					if (P.x == 0.0 && Q.x == 0.0) HalfEdge(3,idx_edge) = -1;  // ghost twin for x=0 face
					if (P.x == 1.0 && Q.x == 1.0) HalfEdge(3,idx_edge) = -2;  // ghost twin for x=1 face
					if (P.y == 0.0 && Q.y == 0.0) HalfEdge(3,idx_edge) = -3;  // ghost twin for y=0 face
					if (P.y == 1.0 && Q.y == 1.0) HalfEdge(3,idx_edge) = -4;  // ghost twin for y=1 face
					if (P.z == 0.0 && Q.z == 0.0) HalfEdge(3,idx_edge) = -5;  // ghost twin for z=0 face
					if (P.z == 1.0 && Q.z == 1.0) HalfEdge(3,idx_edge) = -6;  // ghost twin for z=1 face
				}
				// Find all the angles
				for (int idx=0; idx<EdgeCount; idx++){
					int V1=HalfEdge(0,idx);
					int V2=HalfEdge(1,idx);
					int T1= HalfEdge(2,idx_edge);
					int twin=HalfEdge(3,idx_edge); 
					if (twin == -1){
						
					}
				}

				// Map vertices to global coordinates
				for (int idx=0;idx<VertexCount;idx++) {
					P = cellvertices(idx);
					P.x  += i;
					P.y  += j;
					P.z  += k;
					cellvertices(idx) = P;
				}
			}
		}
	}

}




/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University
  Copyright Equnior ASA

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
#include "analysis/dcel.h"

DCEL::DCEL() {}

DCEL::~DCEL() {
    TriangleCount = 0;
    VertexCount = 0;
}

int DCEL::Face(int index) { return FaceData[index]; }

void DCEL::Write() {
    int e1, e2, e3;
    FILE *TRIANGLES;
    TRIANGLES = fopen("triangles.stl", "w");
    fprintf(TRIANGLES, "solid \n");
    for (int idx = 0; idx < TriangleCount; idx++) {
        e1 = Face(idx);
        e2 = halfedge.next(e1);
        e3 = halfedge.next(e2);
        auto P1 = vertex.coords(halfedge.v1(e1));
        auto P2 = vertex.coords(halfedge.v1(e2));
        auto P3 = vertex.coords(halfedge.v1(e3));
        fprintf(TRIANGLES, "vertex %f %f %f\n", P1.x, P1.y, P1.z);
        fprintf(TRIANGLES, "vertex %f %f %f\n", P2.x, P2.y, P2.z);
        fprintf(TRIANGLES, "vertex %f %f %f\n", P3.x, P3.y, P3.z);
    }
    fclose(TRIANGLES);
}

void DCEL::LocalIsosurface(const DoubleArray &A, double value, const int i,
                           const int j, const int k) {
    Point P, Q;
    Point PlaceHolder;
    Point C0, C1, C2, C3, C4, C5, C6, C7;

    Point VertexList[12];
    Point NewVertexList[12];
    int LocalRemap[12];

    Point cellvertices[20];
    std::array<std::array<int, 3>, 20> Triangles;

    // Values from array 'A' at the cube corners
    double CubeValues[8];

    // Points corresponding to cube corners
    C0.x = 0.0;
    C0.y = 0.0;
    C0.z = 0.0;
    C1.x = 1.0;
    C1.y = 0.0;
    C1.z = 0.0;
    C2.x = 1.0;
    C2.y = 1.0;
    C2.z = 0.0;
    C3.x = 0.0;
    C3.y = 1.0;
    C3.z = 0.0;
    C4.x = 0.0;
    C4.y = 0.0;
    C4.z = 1.0;
    C5.x = 1.0;
    C5.y = 0.0;
    C5.z = 1.0;
    C6.x = 1.0;
    C6.y = 1.0;
    C6.z = 1.0;
    C7.x = 0.0;
    C7.y = 1.0;
    C7.z = 1.0;

    CubeValues[0] = A(i, j, k) - value;
    CubeValues[1] = A(i + 1, j, k) - value;
    CubeValues[2] = A(i + 1, j + 1, k) - value;
    CubeValues[3] = A(i, j + 1, k) - value;
    CubeValues[4] = A(i, j, k + 1) - value;
    CubeValues[5] = A(i + 1, j, k + 1) - value;
    CubeValues[6] = A(i + 1, j + 1, k + 1) - value;
    CubeValues[7] = A(i, j + 1, k + 1) - value;
    //printf("Set cube values: %i, %i, %i \n",i,j,k);

    //Determine the index into the edge table which
    //tells us which vertices are inside of the surface
    int CubeIndex = 0;
    if (CubeValues[0] < 0.0f)
        CubeIndex |= 1;
    if (CubeValues[1] < 0.0f)
        CubeIndex |= 2;
    if (CubeValues[2] < 0.0f)
        CubeIndex |= 4;
    if (CubeValues[3] < 0.0f)
        CubeIndex |= 8;
    if (CubeValues[4] < 0.0f)
        CubeIndex |= 16;
    if (CubeValues[5] < 0.0f)
        CubeIndex |= 32;
    if (CubeValues[6] < 0.0f)
        CubeIndex |= 64;
    if (CubeValues[7] < 0.0f)
        CubeIndex |= 128;

    //Find the vertices where the surface intersects the cube
    if (edgeTable[CubeIndex] & 1) {
        P = VertexInterp(C0, C1, CubeValues[0], CubeValues[1]);
        VertexList[0] = P;
        Q = C0;
    }
    if (edgeTable[CubeIndex] & 2) {
        P = VertexInterp(C1, C2, CubeValues[1], CubeValues[2]);
        VertexList[1] = P;
        Q = C1;
    }
    if (edgeTable[CubeIndex] & 4) {
        P = VertexInterp(C2, C3, CubeValues[2], CubeValues[3]);
        VertexList[2] = P;
        Q = C2;
    }
    if (edgeTable[CubeIndex] & 8) {
        P = VertexInterp(C3, C0, CubeValues[3], CubeValues[0]);
        VertexList[3] = P;
        Q = C3;
    }
    if (edgeTable[CubeIndex] & 16) {
        P = VertexInterp(C4, C5, CubeValues[4], CubeValues[5]);
        VertexList[4] = P;
        Q = C4;
    }
    if (edgeTable[CubeIndex] & 32) {
        P = VertexInterp(C5, C6, CubeValues[5], CubeValues[6]);
        VertexList[5] = P;
        Q = C5;
    }
    if (edgeTable[CubeIndex] & 64) {
        P = VertexInterp(C6, C7, CubeValues[6], CubeValues[7]);
        VertexList[6] = P;
        Q = C6;
    }
    if (edgeTable[CubeIndex] & 128) {
        P = VertexInterp(C7, C4, CubeValues[7], CubeValues[4]);
        VertexList[7] = P;
        Q = C7;
    }
    if (edgeTable[CubeIndex] & 256) {
        P = VertexInterp(C0, C4, CubeValues[0], CubeValues[4]);
        VertexList[8] = P;
        Q = C0;
    }
    if (edgeTable[CubeIndex] & 512) {
        P = VertexInterp(C1, C5, CubeValues[1], CubeValues[5]);
        VertexList[9] = P;
        Q = C1;
    }
    if (edgeTable[CubeIndex] & 1024) {
        P = VertexInterp(C2, C6, CubeValues[2], CubeValues[6]);
        VertexList[10] = P;
        Q = C2;
    }
    if (edgeTable[CubeIndex] & 2048) {
        P = VertexInterp(C3, C7, CubeValues[3], CubeValues[7]);
        VertexList[11] = P;
        Q = C3;
    }

    VertexCount = 0;
    for (int idx = 0; idx < 12; idx++)
        LocalRemap[idx] = -1;

    for (int idx = 0; triTable[CubeIndex][idx] != -1; idx++) {
        if (LocalRemap[triTable[CubeIndex][idx]] == -1) {
            NewVertexList[VertexCount] = VertexList[triTable[CubeIndex][idx]];
            LocalRemap[triTable[CubeIndex][idx]] = VertexCount;
            VertexCount++;
        }
    }

    //printf("Found %i vertices \n",VertexCount);

    for (int idx = 0; idx < VertexCount; idx++) {
        P = NewVertexList[idx];
        //P.x  += i;
        //P.y  += j;
        //P.z  += k;
        cellvertices[idx] = P;
    }

    TriangleCount = 0;
    for (int idx = 0; triTable[CubeIndex][idx] != -1; idx += 3) {
        Triangles[TriangleCount][0] = LocalRemap[triTable[CubeIndex][idx + 0]];
        Triangles[TriangleCount][1] = LocalRemap[triTable[CubeIndex][idx + 1]];
        Triangles[TriangleCount][2] = LocalRemap[triTable[CubeIndex][idx + 2]];
        TriangleCount++;
    }
    int nTris = TriangleCount;

    // Now add the local values to the DCEL data structure
    if (nTris > 0) {
        FaceData.resize(TriangleCount);
        //printf("Construct halfedge structure... \n");
        //printf("   Construct %i triangles \n",nTris);
        halfedge.resize(nTris * 3);
        int idx_edge = 0;
        for (int idx = 0; idx < TriangleCount; idx++) {
            int V1 = Triangles[idx][0];
            int V2 = Triangles[idx][1];
            int V3 = Triangles[idx][2];
            FaceData[idx] = idx_edge;
            // first edge: V1->V2
            halfedge.data(0, idx_edge) = V1;           // first vertex
            halfedge.data(1, idx_edge) = V2;           // second vertex
            halfedge.data(2, idx_edge) = idx;          // triangle
            halfedge.data(3, idx_edge) = -1;           // twin
            halfedge.data(4, idx_edge) = idx_edge + 2; // previous edge
            halfedge.data(5, idx_edge) = idx_edge + 1; // next edge
            idx_edge++;
            // second edge: V2->V3
            halfedge.data(0, idx_edge) = V2;           // first vertex
            halfedge.data(1, idx_edge) = V3;           // second vertex
            halfedge.data(2, idx_edge) = idx;          // triangle
            halfedge.data(3, idx_edge) = -1;           // twin
            halfedge.data(4, idx_edge) = idx_edge - 1; // previous edge
            halfedge.data(5, idx_edge) = idx_edge + 1; // next edge
            idx_edge++;
            // third edge: V3->V1
            halfedge.data(0, idx_edge) = V3;           // first vertex
            halfedge.data(1, idx_edge) = V1;           // second vertex
            halfedge.data(2, idx_edge) = idx;          // triangle
            halfedge.data(3, idx_edge) = -1;           // twin
            halfedge.data(4, idx_edge) = idx_edge - 1; // previous edge
            halfedge.data(5, idx_edge) = idx_edge - 2; // next edge
            idx_edge++;
            //printf("   ***tri %i ***edge %i *** \n",idx, idx_edge);
        }
        //printf("  parsing halfedge structure\n");
        int EdgeCount = idx_edge;
        for (int idx = 0; idx < EdgeCount; idx++) {
            int V1 = halfedge.data(0, idx);
            int V2 = halfedge.data(1, idx);
            // Find all the twins within the cube
            for (int jdx = 0; jdx < EdgeCount; jdx++) {
                if (halfedge.data(1, jdx) == V1 &&
                    halfedge.data(0, jdx) == V2) {
                    // this is the pair
                    halfedge.data(3, idx) = jdx;
                    halfedge.data(3, jdx) = idx;
                }
                if (halfedge.data(1, jdx) == V2 &&
                    halfedge.data(0, jdx) == V1 && !(idx == jdx)) {
                    std::printf(
                        "WARNING: half edges with identical orientation! \n");
                }
            }
            // Use "ghost" twins if edge is on a cube face
            P = cellvertices[V1];
            Q = cellvertices[V2];
            if (P.x == 0.0 && Q.x == 0.0)
                halfedge.data(3, idx) = -1; // ghost twin for x=0 face
            if (P.x == 1.0 && Q.x == 1.0)
                halfedge.data(3, idx) = -4; // ghost twin for x=1 face
            if (P.y == 0.0 && Q.y == 0.0)
                halfedge.data(3, idx) = -2; // ghost twin for y=0 face
            if (P.y == 1.0 && Q.y == 1.0)
                halfedge.data(3, idx) = -5; // ghost twin for y=1 face
            if (P.z == 0.0 && Q.z == 0.0)
                halfedge.data(3, idx) = -3; // ghost twin for z=0 face
            if (P.z == 1.0 && Q.z == 1.0)
                halfedge.data(3, idx) = -6; // ghost twin for z=1 face
        }
    }

    // Map vertices to global coordinates
    for (int idx = 0; idx < VertexCount; idx++) {
        P = cellvertices[idx];
        P.x += i;
        P.y += j;
        P.z += k;
        vertex.assign(idx, P);
    }
}

Point DCEL::TriNormal(int edge) {
    Point P, Q, R;
    Point U, V, W;
    double nx, ny, nz, len;
    // at cube faces define outward normal to cube
    if (edge == -1) {
        W.x = -1.0;
        W.y = 0.0;
        W.z = 0.0; // x cube face
    } else if (edge == -2) {
        W.x = 0.0;
        W.y = -1.0;
        W.z = 0.0; // y cube face
    } else if (edge == -3) {
        W.x = 0.0;
        W.y = 0.0;
        W.z = -1.0; // z cube face
    } else if (edge == -4) {
        W.x = 1.0;
        W.y = 0.0;
        W.z = 0.0; // x cube face
    } else if (edge == -5) {
        W.x = 0.0;
        W.y = 1.0;
        W.z = 0.0; // y cube face
    } else if (edge == -6) {
        W.x = 0.0;
        W.y = 0.0;
        W.z = 1.0; // z cube face
    } else {
        // vertices for triange
        int e2 = halfedge.next(edge);
        int e3 = halfedge.next(e2);
        P = vertex.coords(halfedge.v1(edge));
        Q = vertex.coords(halfedge.v1(e2));
        R = vertex.coords(halfedge.v1(e3));
        // edge vectors
        U = Q - P;
        V = R - Q;
        // normal vector
        nx = U.y * V.z - U.z * V.y;
        ny = U.z * V.x - U.x * V.z;
        nz = U.x * V.y - U.y * V.x;
        len = sqrt(nx * nx + ny * ny + nz * nz);
        W.x = nx / len;
        W.y = ny / len;
        W.z = nz / len;
    }
    return W;
}

double DCEL::EdgeAngle(int edge) {
    double angle;
    double dotprod;
    Point P, Q, R; // triangle vertices
    Point U, V, W; // normal vectors
    int e2 = halfedge.next(edge);
    int e3 = halfedge.next(e2);
    P = vertex.coords(halfedge.v1(edge));
    Q = vertex.coords(halfedge.v1(e2));
    R = vertex.coords(halfedge.v1(e3));
    U = TriNormal(edge);
    V = TriNormal(halfedge.twin(edge));
    if (halfedge.twin(edge) < 0) {
        // compute edge normal in plane of cube face
        W = P - Q; // edge tangent vector
        double length = sqrt(W.x * W.x + W.y * W.y + W.z * W.z);
        W.x /= length;
        W.y /= length;
        W.z /= length;
        // edge normal within the plane of the cube face
        double nx = W.y * V.z - W.z * V.y;
        double ny = W.z * V.x - W.x * V.z;
        double nz = W.x * V.y - W.y * V.x;
        length = sqrt(nx * nx + ny * ny + nz * nz);
        // new value for V is this normal vector
        V.x = nx / length;
        V.y = ny / length;
        V.z = nz / length;
        dotprod = U.x * V.x + U.y * V.y + U.z * V.z;
        if (dotprod < 0.f) {
            //printf("negative dot product on face\n");
            dotprod = -dotprod;
            V.x = -V.x;
            V.y = -V.y;
            V.z = -V.z;
        }

        if (dotprod > 1.f)
            dotprod = 1.f;
        if (dotprod < -1.f)
            dotprod = -1.f;
        angle = acos(dotprod);
        /* project onto plane of cube face also works
		W = U - dotprod*V;
		length = sqrt(W.x*W.x+W.y*W.y+W.z*W.z); // for normalization
		dotprod = (U.x*W.x + U.y*W.y + U.z*W.z)/length;
		if (dotprod > 1.f) dotprod=1.f;
		if (dotprod < -1.f) dotprod=-1.f;
		angle = acos(dotprod);
		 */
    } else {
        dotprod = U.x * V.x + U.y * V.y + U.z * V.z;
        if (dotprod > 1.f)
            dotprod = 1.f;
        if (dotprod < -1.f)
            dotprod = -1.f;
        angle = 0.5 * acos(dotprod);
    }
    // determine if angle is concave or convex based on edge normal
    W.x = (P.y - Q.y) * U.z - (P.z - Q.z) * U.y;
    W.y = (P.z - Q.z) * U.x - (P.x - Q.x) * U.z;
    W.z = (P.x - Q.x) * U.y - (P.y - Q.y) * U.x;
    //length = sqrt(nx*nx+ny*ny+nz*nz);
    Point w = 0.5 * (P + Q) - R;
    if (W.x * w.x + W.y * w.y + W.z * w.z < 0.f) {
        //printf("flip edge normal \n");
        W.x = -W.x;
        W.y = -W.y;
        W.z = -W.z;
    }
    if (W.x * V.x + W.y * V.y + W.z * V.z > 0.f) {
        // concave
        angle = -angle;
    }
    if (angle != angle)
        angle = 0.0;
    //printf("angle=%f,dot=%f (Edge=%i, twin=%i): P={%f, %f, %f}, Q={%f, %f, %f} U={%f, %f, %f}, V={%f, %f, %f}\n",angle,dotprod,edge,halfedge.twin(edge),P.x,P.y,P.z,Q.x,Q.y,Q.z,U.x,U.y,U.z,V.x,V.y,V.z);
    return angle;
}

void iso_surface(const Array<double> &Field, const double isovalue) {
    DCEL object;
    int e1, e2, e3;
    FILE *TRIANGLES;
    TRIANGLES = fopen("isosurface.stl", "w");
    fprintf(TRIANGLES, "solid isosurface\n");
    int Nx = Field.size(0);
    int Ny = Field.size(1);
    int Nz = Field.size(2);
    for (int k = 1; k < Nz - 1; k++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                object.LocalIsosurface(Field, isovalue, i, j, k);
                for (int idx = 0; idx < object.TriangleCount; idx++) {
                    e1 = object.Face(idx);
                    e2 = object.halfedge.next(e1);
                    e3 = object.halfedge.next(e2);
                    auto P1 = object.vertex.coords(object.halfedge.v1(e1));
                    auto P2 = object.vertex.coords(object.halfedge.v1(e2));
                    auto P3 = object.vertex.coords(object.halfedge.v1(e3));
                    auto Normal = object.TriNormal(e1);
                    // P1.x += 1.0*i; P1.y += 1.0*j; P1.z +=1.0*k;
                    //P2.x += 1.0*i; P2.y += 1.0*j; P2.z +=1.0*k;
                    //P3.x += 1.0*i; P3.y += 1.0*j; P3.z +=1.0*k;
                    fprintf(TRIANGLES, "facet normal %f %f %f\n", Normal.x,
                            Normal.y, Normal.z);
                    fprintf(TRIANGLES, "   outer loop\n");
                    fprintf(TRIANGLES, "      vertex %f %f %f\n", P1.x, P1.y,
                            P1.z);
                    fprintf(TRIANGLES, "      vertex %f %f %f\n", P2.x, P2.y,
                            P2.z);
                    fprintf(TRIANGLES, "      vertex %f %f %f\n", P3.x, P3.y,
                            P3.z);
                    fprintf(TRIANGLES, "   endloop\n");
                    fprintf(TRIANGLES, "endfacet\n");
                }
            }
        }
    }
    fprintf(TRIANGLES, "endsolid isosurface\n");
    fclose(TRIANGLES);
}

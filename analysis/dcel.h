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
#ifndef DCEL_INC
#define DCEL_INC

#include <vector>
#include "analysis/pmmc.h"

/**
 * \class Vertex
 * @brief store vertex for DCEL data structure
*/

// Vertex structure
class Vertex {
public:
    Vertex() { d_data.resize(12); }
    ~Vertex() = default;
    Vertex(const Vertex &) = delete;
    Vertex operator=(const Vertex &) = delete;

    // Add/assign a point
    inline void add(const Point &P) { d_data.push_back(P); }
    inline void assign(int idx, const Point &P) { d_data[idx] = P; }

    // Get a point
    inline Point &coords(int idx) { return d_data[idx]; }
    inline const Point &coords(int idx) const { return d_data[idx]; }

    int IncidentEdge();

    // Return the number of points
    inline int size() const { return d_data.size(); }

private:
    std::vector<Point> d_data;
};

/**
 * \class Halfedge
 * @brief store half edge for DCEL data structure
*/
class Halfedge {
public:
    Halfedge() = default;
    ~Halfedge() = default;
    Halfedge(const Halfedge &) = delete;
    Halfedge operator=(const Halfedge &) = delete;

    inline int v1(int edge) const { return d_data[edge][0]; }
    inline int v2(int edge) const { return d_data[edge][1]; }
    inline int face(int edge) const { return d_data[edge][2]; }
    inline int twin(int edge) const { return d_data[edge][3]; }
    inline int prev(int edge) const { return d_data[edge][4]; }
    inline int next(int edge) const { return d_data[edge][5]; }

    inline int size() const { return d_data.size(); }
    inline void resize(int N) { d_data.resize(N); }

    inline int &data(int i, int j) { return d_data[j][i]; }
    inline const int &data(int i, int j) const { return d_data[j][i]; }

private:
    std::vector<std::array<int, 6>> d_data;
};

/**
 * \class DCEL
 * @details doubly connected edge list data structure 
*/
class DCEL {
public:
    DCEL();
    ~DCEL();

    int face();
    Vertex vertex;
    Halfedge halfedge;
    void LocalIsosurface(const DoubleArray &A, double value, int i, int j,
                         int k);
    void Write();
    int Face(int index);

    double origin(int edge);
    double EdgeAngle(int edge);
    Point TriNormal(int edge);
    int TriangleCount;
    int VertexCount;

private:
    std::vector<int> FaceData;
};

void iso_surface(const Array<double> &Field, const double isovalue);

#endif

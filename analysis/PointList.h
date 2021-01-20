/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University

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
#ifndef PointList_INC
#define PointList_INC

#include <math.h>

struct LBPM_Point {
    LBPM_Point() : x(0.0), y(0.0), z(0.0) {}
    LBPM_Point(double xv,double yv,double zv) : x(xv), y(yv), z(zv) {}
    LBPM_Point(const LBPM_Point& rhs): x(rhs.x), y(rhs.y), z(rhs.z) {}
    //Point& operator=(const Point& rhs) { this->x=rhs.x; this->y=rhs.y; this->z=rhs.z; return *this; }
    //~Point() {}
    double x,y,z;
};
typedef LBPM_Point Point;

inline Point operator+(const Point &A,const Point &B) {return Point(A.x+B.x,A.y+B.y,A.z+B.z);}
inline Point operator-(const Point &A,const Point &B) {return Point(A.x-B.x,A.y-B.y,A.z-B.z);}
inline Point operator*(const Point &A,double v) {return Point(A.x*v,A.y*v,A.z*v);}
inline Point operator*(double v,const Point &A) {return Point(A.x*v,A.y*v,A.z*v);}
inline Point operator/(const Point &A,double v) {return Point(A.x/v,A.y/v,A.z/v);}
inline Point operator-(const Point &A) {return Point(-A.x,-A.y,-A.z);}

inline bool operator==(const Point &A,const Point &B) {return (A.x==B.x && A.y==B.y && A.z==B.z);}
inline bool operator!=(const Point &A,const Point &B) {return (A.x!=B.x || A.y!=B.y || A.z!=B.z);}

inline double Norm(const Point &A) {return sqrt(A.x*A.x+A.y*A.y+A.z*A.z);}
inline Point Cross(const Point &A,const Point &B) {return Point(A.y*B.z-A.z*B.y,B.x*A.z-A.x*B.z,A.x*B.y-A.y*B.x);}
inline double Dot(const Point &A,const Point &B) {return (A.x*B.x+A.y*B.y+A.z*B.z);}
inline double Distance(const Point &A,const Point &B) {return sqrt((A.x-B.x)*(A.x-B.x)+(A.y-B.y)*(A.y-B.y)+(A.z-B.z)*(A.z-B.z));}

/*
class  PointList{
public:
	int length;
	int m;
	int index;

	Point pt;

	double *data;
	PointList();
	PointList(int size);
	~PointList();

	void New(int size);

	Point & operator()(int idx)
	{
		index = idx;
		pt.x = data[3*index];
		pt.y = data[3*index+1];
		pt.z = data[3*index+2];
		return pt;
	}
	
	Point &operator=(Point &P){
		pt.x = P.x;
		pt.y = P.y;
		pt.z = P.z;
		data[3*index]=pt.x;
		data[3*index+1]=pt.y;
		data[3*index+2]=pt.z;
	}
};

// *****************************************
// ******** class PointList **************
// *****************************************
PointList::PointList()
{
	m=length=0;
}

PointList::PointList(int size)
{
	m=size;
	data = new double [3*size];
	length = size;
}

void PointList::New(int size)
{
	m=size;
	data = new double [3*size];
	length = size;
}


PointList::~PointList()
{
	delete data;
}
*/
template <class T>
class DTList {
public:
    DTList() : Data(0), length(0), refCount(new size_t(1)), outOfRange() {}
    DTList(const DTList<T> &A) : Data(A.Data), length(A.length), refCount(A.refCount), outOfRange() {++(*refCount);}
protected:
    DTList(size_t len) : Data(len<=0 ? 0 : new T[len]), length(len<=0 ? 0 : len), refCount(new size_t(1)), outOfRange() {}
public:
    
    virtual ~DTList() {
      --(*refCount);
      if (*refCount==0) {delete [] Data; delete refCount;}
      Data = 0; refCount = 0; length=0;
    }
	
    DTList<T> &operator=(const DTList<T> &A) {
        if (A.refCount!=refCount) { // Otherwise doing A=A.
            --(*refCount);
            if (*refCount==0) {delete [] Data; delete refCount;}
            refCount = A.refCount;
            ++(*refCount);
            length = A.length;
            Data = A.Data;
        }
        return *this;
    }
	
    size_t MemoryUsed(void) const {return length*sizeof(T);}
	
    const T *Pointer(void) const {return Data;}
	
    size_t IsEmpty(void) const {return (Data==0);}
    size_t Length(void) const {return length;}

    const T operator()(size_t i) const  {return Data[i];}
	
protected:
    T *Data;
    size_t length;
    size_t *refCount;
	
    // Should be static.
    T outOfRange;
};

template <class T>
class DTMutableList : public DTList<T> {
public:
    DTMutableList() : DTList<T>() {}
    DTMutableList(size_t len) : DTList<T>(len) {}
    DTMutableList(const DTMutableList<T> &A) : DTList<T>(A) {}
	
    DTMutableList<T> &operator=(const DTMutableList<T> &A) {DTList<T>::operator=(A); return *this;}
	
    T *Pointer(void) {return DTList<T>::Data;}
    const T *Pointer(void) const {return DTList<T>::Data;}
    T &operator()(size_t i) {return DTList<T>::Data[i];}
    T operator()(size_t i) const  {return DTList<T>::Data[i];}
	
    DTMutableList<T> &operator=(T v) {for (size_t i=0;i<DTList<T>::length;i++) DTList<T>::Data[i] = v; return *this;}
};

template <class T> DTMutableList<T> TruncateSize(const DTList<T> &A,size_t length)
{
    if (length>A.Length()) length = A.Length();
    DTMutableList<T> toReturn(length);
    const T *fromP = A.Pointer();
    T *toP = toReturn.Pointer();
    for (size_t i=0;i<length;i++) toP[i] = fromP[i];
    return toReturn;
}

template <class T> DTMutableList<T> IncreaseSize(const DTList<T> &A,size_t addLength)
{
    DTMutableList<T> toReturn(A.Length()+(addLength>=0 ? addLength : 0));
    size_t len = A.Length();
    const T *fromP = A.Pointer();
    T *toP = toReturn.Pointer();
    for (size_t i=0;i<len;i++) toP[i] = fromP[i];
    return toReturn;
}

#endif



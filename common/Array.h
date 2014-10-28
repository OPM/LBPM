#ifndef ARRAY_H_INC
#define ARRAY_H_INC

#include <iostream>
#include <string.h>

// ********** ARRAY CLASS INFO **************************************
/*
 //..............................................................
 // Overview of array classes in Array.h
 //......... Declaration of array objects........................
 // supports integer data (IntArray) or double data (DoubleArray)
 // supports one, two or three dimensional arrays
 IntArray A(m);			// array of integers, Length n
 IntArray A(m,n);		// array of integers, dimensions: m x n
 DoubleArray A(m,n,o);	// array of doubles, dimensions m x n x o	
 //............ Access the size of the array.....................
 A.m;				// size of first dimension
 A.n;				// size of second dimension
 A.o;				// size of third dimension
 A.Length;			// total number of values stored
 //........... Access the array entries .........................
 A(i);				// return data[i] 
 A(i,j);			// return data[j*m+i]
 A(i,j,k);			// return data[k*m*n+j*m+i]
 //..............................................................
 */

using namespace std;

class  IntArray{
public:
	IntArray Copy();

	int Length;
	int m,n,o;
	int *data;
	IntArray();
	IntArray(int size);
	IntArray(int nx,int ny);
	IntArray(int nx,int ny,int nz);
	~IntArray();
	
	void New(int size);
	void New(int nx, int ny);
	void New(int nx, int ny, int nz);
	
	int & operator()(int index)
	{return data[index];}	
	int & operator()(int i, int j)
	{ return data[j*m+i];}
	int & operator()(int i, int j, int k)
	{ return data[k*m*n+j*m+i];}
	
	int e(int i);
	int e(int i, int j);
	int e(int i, int j, int k);
	
	int *Pointer() {return data;}
};

class  DoubleArray{
public:
	DoubleArray Copy();

	int Length;
	int m;
	int n;
	int o;
	double *data;
	DoubleArray();
	DoubleArray(int size);
	DoubleArray(int nx, int ny);
	DoubleArray(int nx, int ny, int nz);
	~DoubleArray();
	
	void New(int size);
	void New(int nx, int ny);
	void New(int nx, int ny, int nz);
	
	double & operator()(int index)
	{return data[index];}	
	double & operator()(int i, int j)
	{return data[j*m+i];}
	double & operator()(int i, int j, int k)
	{return data[k*m*n+j*m+i];}
	double e(int i);
	double e(int i, int j);
	double e(int i, int j, int k);

	double *Pointer() {return data;}
};

/*
 *  Array.cpp
 *
 *  Created by James Mcclure on 3/31/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

// *****************************************
// ******** class IntArray  ****************
// *****************************************
/*IntArray::IntArray();
IntArray::IntArray(int size);
IntArray::IntArray(int nx, int ny);
IntArray::IntArray(int nx, int ny, int nz);
IntArray::~IntArray();
void IntArray::New(int size);
void IntArray::New(int nx, int ny);
void IntArray::New(int nx, int ny,int nz);
int IntArray::e(int i);
int IntArray::e(int i, int j);
int IntArray::e(int i, int j, int k);

// *****************************************
// ******** class DoubleArray **************
// *****************************************
DoubleArray::DoubleArray();
DoubleArray::DoubleArray(int size);
DoubleArray::DoubleArray(int nx, int ny);
DoubleArray::DoubleArray(int nx, int ny,int nz);
DoubleArray::~DoubleArray();
void DoubleArray::New(int size);
void DoubleArray::New(int nx, int ny);
void DoubleArray::New(int nx, int ny,int nz);
double DoubleArray::e(int i);
double DoubleArray::e(int i, int j);
double DoubleArray::e(int i, int j, int k);
extern DoubleArray IncreaseSize(DoubleArray &A, int addLength);
extern IntArray IncreaseSize(IntArray &A, int addLength);
*/

/*
 *  Array.cpp
 *
 *  Created by James Mcclure on 3/31/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
/*
 *  Array.cpp
 *
 *  Created by James Mcclure on 3/31/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

IntArray::IntArray()
{
	m=n=o=Length=0;
}

IntArray::IntArray(int size)
{
	data = new int [size];
	Length = size;
	m = size;
	n = 1;
	o = 1;
}

IntArray::IntArray(int nx, int ny)
{
	m = nx;
	n = ny;
	o = 1;
	Length = m*n;
	data = new int [Length];
}

IntArray::IntArray(int nx, int ny, int nz)
{
	m = nx;
	n = ny;
	o = nz;
	Length = m*n*o;
	data = new int [Length];
}

IntArray::~IntArray()
{
	delete data;
}


void IntArray::New(int size)
{
	m=size;
	n = 1;
	o = 1;
	data = new int [size];
	Length = size;
}

void IntArray::New(int nx, int ny)
{
	m = nx;
	n = ny;
	o = 1;
	Length = m*n;
	data = new int [Length];
}

void IntArray::New(int nx, int ny,int nz)
{
	m = nx;
	n = ny;
	o = nz;
	Length = m*n*o;
	data = new int [Length];
}

int IntArray::e(int i)
{
	return data[i];
}

int IntArray::e(int i, int j)
{
	return data[m*j+i];
}

int IntArray::e(int i, int j, int k)
{
	return data[m*n*k+m*j+i];
}

// *****************************************
// ******** class DoubleArray **************
// *****************************************
DoubleArray::DoubleArray()
{
	m=n=o=Length=0;
}

DoubleArray::DoubleArray(int size)
{
	m=size;
	n = 1;
	o = 1;
	data = new double [size];
	Length = size;
}

DoubleArray::DoubleArray(int nx, int ny)
{
	m = nx;
	n = ny;
	o = 1;
	Length = m*n;
	data = new double [Length];
}

DoubleArray::DoubleArray(int nx, int ny,int nz)
{
	m = nx;
	n = ny;
	o = nz;
	Length = m*n*o;
	data = new double [Length];
}

void DoubleArray::New(int size)
{
	m=size;
	n = 1;
	o = 1;
	data = new double [size];
	Length = size;
}

void DoubleArray::New(int nx, int ny)
{
	m = nx;
	n = ny;
	o = 1;
	Length = m*n;
	data = new double [Length];
}

void DoubleArray::New(int nx, int ny,int nz)
{
	m = nx;
	n = ny;
	o = nz;
	Length = m*n*o;
	data = new double [Length];
}

DoubleArray::~DoubleArray()
{
	delete data;
}

double DoubleArray::e(int i)
{
	return data[i];
}	

double DoubleArray::e(int i, int j)
{
	return data[j*m+i];
}

double DoubleArray::e(int i, int j, int k)
{
	return data[k*m*n+j*m+i];
}

extern DoubleArray IncreaseSize(DoubleArray &A, int addLength)
{
    if (addLength<0) {
        printf("IncreaseSize(Array,Length) Length needs to be >0.");
        return DoubleArray();
    }

    int newM,newN,newO;
    if (A.o>1) {
        if (addLength%(A.m*A.n)!=0) {
        	printf("IncreaseSize(Array,Length) Length needs to be a multiple of m*n");
            return DoubleArray();
        }
        newM = A.m;
        newN = A.n;
        newO = A.o + addLength/(A.m*A.n);
    }
    else if (A.n>1) {
        if (addLength%(A.m)!=0) {
        	printf("IncreaseSize(Array,Length) Length needs to be a multiple of m");
            return DoubleArray();
        }
        newM = A.m;
        newN = A.n + addLength/A.m;
        newO = 1;
    }
    else {
        newM = A.m + addLength;
        newN = 1;
        newO = 1;
    }

    DoubleArray toReturn(newM,newN,newO);
    memcpy(toReturn.Pointer(),A.Pointer(),A.Length*sizeof(double));
    return toReturn;
}
extern IntArray IncreaseSize(IntArray &A, int addLength)
{
    if (addLength<0) {
        printf("IncreaseSize(Array,Length) Length needs to be >0.");
        return IntArray();
    }

    int newM,newN,newO;
    if (A.o>1) {
        if (addLength%(A.m*A.n)!=0) {
        	printf("IncreaseSize(Array,Length) Length needs to be a multiple of m*n");
            return IntArray();
        }
        newM = A.m;
        newN = A.n;
        newO = A.o + addLength/(A.m*A.n);
    }
    else if (A.n>1) {
        if (addLength%(A.m)!=0) {
        	printf("IncreaseSize(Array,Length) Length needs to be a multiple of m");
            return IntArray();
        }
        newM = A.m;
        newN = A.n + addLength/A.m;
        newO = 1;
    }
    else {
        newM = A.m + addLength;
        newN = 1;
        newO = 1;
    }

    IntArray toReturn(newM,newN,newO);
    memcpy(toReturn.Pointer(),A.Pointer(),A.Length*sizeof(int));
    return toReturn;
}

DoubleArray DoubleArray::Copy()
{
    DoubleArray CopyInto(m,n,o);
    // Check that the allocation worked.
    if (CopyInto.Length!=Length) return CopyInto; // Failed.  Already printed an error message.
    memcpy(CopyInto.Pointer(),Pointer(),Length*sizeof(double));
    return CopyInto;
}

#endif

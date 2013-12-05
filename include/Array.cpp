#if 0
#include "Array.h"
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
        printf("IncreaseSize(Array,Length)","Length needs to be >0.");
        return DoubleArray();
    }

    int newM,newN,newO;
    if (A.o>1) {
        if (addLength%(A.m*A.n)!=0) {
        	printf("IncreaseSize(Array,Length)","Length needs to be a multiple of m*n");
            return DoubleArray();
        }
        newM = A.m;
        newN = A.n;
        newO = A.o + addLength/(A.m*A.n);
    }
    else if (A.n>1) {
        if (addLength%(A.m)!=0) {
        	printf("IncreaseSize(Array,Length)","Length needs to be a multiple of m");
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
        printf("IncreaseSize(Array,Length)","Length needs to be >0.");
        return IntArray();
    }

    int newM,newN,newO;
    if (A.o>1) {
        if (addLength%(A.m*A.n)!=0) {
        	printf("IncreaseSize(Array,Length)","Length needs to be a multiple of m*n");
            return IntArray();
        }
        newM = A.m;
        newN = A.n;
        newO = A.o + addLength/(A.m*A.n);
    }
    else if (A.n>1) {
        if (addLength%(A.m)!=0) {
        	printf("IncreaseSize(Array,Length)","Length needs to be a multiple of m");
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

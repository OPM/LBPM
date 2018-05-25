#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#define MAX_LOCAL_BLOB_COUNT 500

using namespace std;

struct BlobInfo{
	
	BlobInfo(int size){
		Values = new double [40*size];
		Ranks = new int [19*size];
		RankLabel = new int [19*size];
	};
	double * Values;
	int * Ranks;
	int * RankLabel;
	// local rank and blob_id

	int rank(int i){return Ranks[19*i];}
	int rank_x(int i){return Ranks[19*i+1];}
	int rank_y(int i){return Ranks[19*i+2];}
	int rank_z(int i){return Ranks[19*i+3];}
	int rank_X(int i){return Ranks[19*i+4];}
	int rank_Y(int i){return Ranks[19*i+5];}
	int rank_Z(int i){return Ranks[19*i+6];}
	int rank_xy(int i){return Ranks[19*i+7];}
	int rank_xY(int i){return Ranks[19*i+8];}
	int rank_Xy(int i){return Ranks[19*i+9];}
	int rank_XY(int i){return Ranks[19*i+10];}
	int rank_xz(int i){return Ranks[19*i+11];}
	int rank_xZ(int i){return Ranks[19*i+12];}
	int rank_Xz(int i){return Ranks[19*i+13];}
	int rank_XZ(int i){return Ranks[19*i+14];}
	int rank_yz(int i){return Ranks[19*i+15];}
	int rank_yZ(int i){return Ranks[19*i+16];}
	int rank_Yz(int i){return Ranks[19*i+17];}
	int rank_YZ(int i){return Ranks[19*i+18];}	
	
	int blobID(int i){return RankLabel[19*i];}
	int blobID_x(int i){return RankLabel[19*i+1];}
	int blobID_y(int i){return RankLabel[19*i+2];}
	int blobID_z(int i){return RankLabel[19*i+3];}
	int blobID_X(int i){return RankLabel[19*i+4];}
	int blobID_Y(int i){return RankLabel[19*i+5];}
	int blobID_Z(int i){return RankLabel[19*i+6];}
	int blobID_xy(int i){return RankLabel[19*i+7];}
	int blobID_xY(int i){return RankLabel[19*i+8];}
	int blobID_Xy(int i){return RankLabel[19*i+9];}
	int blobID_XY(int i){return RankLabel[19*i+10];}
	int blobID_xz(int i){return RankLabel[19*i+11];}
	int blobID_xZ(int i){return RankLabel[19*i+12];}
	int blobID_Xz(int i){return RankLabel[19*i+13];}
	int blobID_XZ(int i){return RankLabel[19*i+14];}
	int blobID_yz(int i){return RankLabel[19*i+15];}
	int blobID_yZ(int i){return RankLabel[19*i+16];}
	int blobID_Yz(int i){return RankLabel[19*i+17];}
	int blobID_YZ(int i){return RankLabel[19*i+18];}	
	
	double voln(int i){return Values[40*i];}
	double pn(int i){return Values[40*i+1];}
	double vnx(int i){return Values[40*i+2];}
	double vny(int i){return Values[40*i+3];}
	double vnz(int i){return Values[40*i+4];}
	double vwnx(int i){return Values[40*i+5];}
	double vwny(int i){return Values[40*i+6];}
	double vwnz(int i){return Values[40*i+7];}
	double Jwn(int i){return Values[40*i+8];}
	double Kwn(int i){return Values[40*i+9];}
	double awn(int i){return Values[40*i+10];}
	double ans(int i){return Values[40*i+11];}
	double Gwnxx(int i){return Values[40*i+12];}
	double Gwnyy(int i){return Values[40*i+13];}
	double Gwnzz(int i){return Values[40*i+14];}
	double Gwnxy(int i){return Values[40*i+15];}
	double Gwnxz(int i){return Values[40*i+16];}
	double Gwnyz(int i){return Values[40*i+17];}
	double Gnsxx(int i){return Values[40*i+18];}
	double Gnsyy(int i){return Values[40*i+19];}
	double Gnszz(int i){return Values[40*i+20];}
	double Gnsxy(int i){return Values[40*i+21];}
	double Gnsxz(int i){return Values[40*i+22];}
	double Gnsyz(int i){return Values[40*i+23];}
	double lwns(int i){return Values[40*i+24];}
	double cospwns(int i){return Values[40*i+25];}
	double vwnsx(int i){return Values[40*i+26];}
	double vwnsy(int i){return Values[40*i+27];}
	double vwnsz(int i){return Values[40*i+28];}
	double cx(int i){return Values[40*i+29];}
	double cy(int i){return Values[40*i+30];}
	double cz(int i){return Values[40*i+31];}

};

int main(void){
	
	int Nx,Ny,Nz,N;
	Nx = Ny = Nz = 100;
	N = Nx*Ny*Nz;
			
	BlobInfo Blobs(MAX_LOCAL_BLOB_COUNT);
//	LocalBlobs = new BlobInfo[MAX_LOCAL_BLOB_COUNT];
	
	double * spheres;
	spheres = new double [N];
	
	int nblobs = 10;
	
	for (int i=0; i<nblobs; i++){
		Blobs.rank(i) = 0;
	}
	
}

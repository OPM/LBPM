#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <Array.h>

using namespace std;
//--------------------------------------------------------------------------------------------------------

inline double PhaseAverageScalar(int *PhaseID, double *Scalar, int N, int Np)
{
	int i;
	// Store the averaged values for each phase
	double *PhaseAveragedValues;
	double *PhaseVolumes;
	PhaseAveragedValues = new double[Np];
	for (i=0; i<N; i++){

	}
}
inline int AssignGlobalBlobID(IntArray &blobs, int &nblobs, int &ncubes, IntArray &indicator, &IntArray cubes_in_blob, int number_of_blobs_below_rank)
{
	// This routine takes the local blob ID and assigns a globally unique tag to each blob
	// the local tag is already unique, we just need to add it to the number of blobs determined
	// for all ranks < the rank of the processor calling this routine
	// Note that the total number of connected blobs is less than number_of_blobs_below_rank
	// since the tag that results after calling this routine does not reflect connectivity 
	// between different MPI sub-domains. 
	int i,j,k;
	int start=0,finish;
	int a,c;

	int minimumBlobID;
	
	for (a=0;a<nblobs;a++){
		finish = start+cubes_in_blob(a);
		//...............
		for (c=start;c<finish;c++){
			// Get cube from the list
			i = blobs(0,c);
			j = blobs(1,c);
			k = blobs(2,c);
			// reassign a unique ID for the whole blob
			indicator(i,j,k) += number_of_blobs_below_rank;
		}
		start = finish;
	}
}

inline int ReassignLocalBlobID(IntArray &blobs, int &nblobs, int &ncubes, IntArray &indicator, &IntArray cubes_in_blob)
{
	// This routine reassigns all blobs within given MPI process so that a unique tag is assigned to
	// all cubes within the blob. If a blob is entirely within the sub-domain, the blob will already 
	// be uniquely tagged. If the blob crosses one or more processor boundaries, multiple tags will be 
	// assigned and this routine will assign the minimum tag found within the blob boundaries.
	// The routine returns the number of blobs within the MPI sub-domain that have had ID reassignments. 
	// To explore the final connectivity, this routine must be run iteratively until no MPI
	// process reassigns the label for a cube.
	int i,j,k;
	int start=0,finish;
	int a,c;
	int count_reassigned_blobs=0;

	int minimumBlobID;
	
	for (a=0;a<nblobs;a++){
		finish = start+cubes_in_blob(a);
		//...............
		i = blobs(0,start);
		j = blobs(1,start);
		k = blobs(2,start);
		//...............
		minimumBlobID = indicator(i,j,k);
		//...............
		// There may be multiple blob IDs for one blob due to parallel assignment
		// Find the minimum blob ID within this blob 
		for (c=start;c<finish;c++){
			// Get cube from the list
			i = blobs(0,c);
			j = blobs(1,c);
			k = blobs(2,c);			
			if (IntArray(i,j,k) < minimumBlobID){
				minimumBlobID = indicator(i,j,k);
				count_reassigned_blobs++;
			}
		}
		for (c=start;c<finish;c++){
			// Get cube from the list
			i = blobs(0,c);
			j = blobs(1,c);
			k = blobs(2,c);
			// reassign a unique ID for the whole blob
			indicator(i,j,k) = minimumBlobID;
		}
		start = finish;
	}
	return count_reassigned_blobs;
}

inline int ComputeLocalBlob(IntArray blobs, int &nblobs, int &ncubes, IntArray indicator,
					   DoubleArray F, DoubleArray S, double vf, double vs, int startx, int starty,
					   int startz, IntArray temp)
{
	// Compute the blob (F>vf|S>vs) starting from (i,j,k) - oil blob
	// F>vf => oil phase S>vs => in porespace
	// update the list of blobs, indicator mesh
	int m = F.m;  // maxima for the meshes
	int n = F.n;
	int o = F.o;

	int cubes_in_blob=0;
	int nrecent = 1;						// number of nodes added at most recent sweep
	temp(0,0) = startx;				// Set the initial point as a "seed" for the sweeps
	temp(1,0) = starty;
	temp(2,0) = startz;
	int ntotal = 1;					// total number of nodes in blob
	indicator(startx,starty,startz) = nblobs;

	int p,s,x,y,z,start,finish,nodx,nody,nodz;
	int imin=startx,imax=startx,jmin=starty,jmax=starty;	// initialize maxima / minima
	int kmin=startz,kmax=startz;
	int d[26][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},
	{1,1,0},{1,-1,0},{-1,1,0},{-1,-1,0},{1,0,1},{-1,0,1},
	{1,0,-1},{-1,0,-1},{0,1,1},{0,-1,1},{0,1,-1},{0,-1,-1},
	{1,1,1},{1,1,-1},{1,-1,1},{1,-1,-1},{-1,1,1},{-1,1,-1},
	{-1,-1,1},{-1,-1,-1}};   // directions to neighbors
	int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};  // cube corners
	bool status = 1;						// status == true => continue to look for points
	while (status == 1){
		start = ntotal - nrecent;
		finish = ntotal;
		nrecent = 0;						// set recent points back to zero for next sweep through
		for (s=start;s<finish;s++){
			// Loop over recent points; look for new points
			x = temp(0,s);
			y = temp(1,s);
			z = temp(2,s);
			// Looop over the directions
			for (p=0;p<26;p++){
				nodx=x+d[p][0];
				if (nodx < 0 ){ nodx = m-1; }		// Periodic BC for x
				if (nodx > m-1 ){ nodx = 0; }
				nody=y+d[p][1];
				if (nody < 0 ){ nody = n-1; }	// Periodic BC for y
				if (nody > n-1 ){ nody = 0; }
				nodz=z+d[p][2];
				if (nodz < 0 ){ nodz = 0; }		// No periodic BC for z
				if (nodz > o-1 ){ nodz = o-1; }
				if ( F(nodx,nody,nodz) > vf && S(nodx,nody,nodz) > vs
					 && indicator(nodx,nody,nodz) == -1 ){
					// Node is a part of the blob - add it to the list
					temp(0,ntotal) = nodx;
					temp(1,ntotal) = nody;
					temp(2,ntotal) = nodz;
					ntotal++;
					nrecent++;
					// Update the indicator map
					indicator(nodx,nody,nodz) = nblobs;
					// Update the min / max for the cube loop
					if ( nodx < imin ){ imin = nodx; }
					if ( nodx > imax ){ imax = nodx; }
					if ( nody < jmin ){ jmin = nody; }
					if ( nody > jmax ){ jmax = nody; }
					if ( nodz < kmin ){ kmin = nodz; }
					if ( nodz > kmax ){ kmax = nodz; }
				}
				else if (F(nodx,nody,nodz) > vf && S(nodx,nody,nodz) > vs
						 && indicator(nodx,nody,nodz) > -1 &&  indicator(nodx,nody,nodz) != nblobs){
					// Some kind of error in algorithm
					printf("Error in blob search algorithm!");
				}
			}

		}
		if ( nrecent == 0){
			status = 0;
		}
	}
	// Use points in temporary storage array to add cubes to the list of blobs
	if ( imin > 0) { imin = imin-1; }
//	if ( imax < m-1) { imax = imax+1; }
	if ( jmin > 0) { jmin = jmin-1; }
//	if ( jmax < n-1) { jmax = jmax+1; }
	if ( kmin > 0) { kmin = kmin-1; }
//	if ( kmax < o-1) { kmax = kmax+1; }
	int i,j,k;
	bool add;
	for (k=kmin;k<kmax;k++){
		for (j=jmin;j<jmax;j++){
			for (i=imin;i<imax;i++){
				// If cube(i,j,k) has any nodes in blob, add it to the list
				// Loop over cube edges
				add = 0;
				for (p=0;p<8;p++){
					nodx = i+cube[p][0];
					nody = j+cube[p][1];
					nodz = k+cube[p][2];
					if ( indicator(nodx,nody,nodz) == nblobs ){
						// Cube corner is in this blob
						add = 1;
					}
				}
				if (add == 1){
					// Add cube to the list
					blobs(0,ncubes) = i;
					blobs(1,ncubes) = j;
					blobs(2,ncubes) = k;
					ncubes++;
					cubes_in_blob++;
					// Loop again to check for overlap
					for (p=0;p<8;p++){
						nodx = i+cube[p][0];
						nody = j+cube[p][1];
						nodz = k+cube[p][2];
						if (indicator(nodx,nody,nodz) > -1 && indicator(nodx,nody,nodz) != nblobs){
							printf("Overlapping cube!");
							cout << i << ", " << j << ", " << k << endl;
						}
					}
				}
			}
		}
	}

	return cubes_in_blob;
}

inline int ComputeBlob(IntArray blobs, int &nblobs, int &ncubes, IntArray indicator,
					   DoubleArray F, DoubleArray S, double vf, double vs, int startx, int starty,
					   int startz, IntArray temp)
{
	// Compute the blob (F>vf|S>vs) starting from (i,j,k) - oil blob
	// F>vf => oil phase S>vs => in porespace
	// update the list of blobs, indicator mesh
	int m = F.m;  // maxima for the meshes
	int n = F.n;
	int o = F.o;

	int cubes_in_blob=0;
	int nrecent = 1;						// number of nodes added at most recent sweep
	temp(0,0) = startx;				// Set the initial point as a "seed" for the sweeps
	temp(1,0) = starty;
	temp(2,0) = startz;
	int ntotal = 1;					// total number of nodes in blob
	indicator(startx,starty,startz) = nblobs;

	int p,s,x,y,z,start,finish,nodx,nody,nodz;
	int imin=startx,imax=startx,jmin=starty,jmax=starty;	// initialize maxima / minima
	int kmin=startz,kmax=startz;
	int d[26][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},
	{1,1,0},{1,-1,0},{-1,1,0},{-1,-1,0},{1,0,1},{-1,0,1},
	{1,0,-1},{-1,0,-1},{0,1,1},{0,-1,1},{0,1,-1},{0,-1,-1},
	{1,1,1},{1,1,-1},{1,-1,1},{1,-1,-1},{-1,1,1},{-1,1,-1},
	{-1,-1,1},{-1,-1,-1}};   // directions to neighbors
	int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};  // cube corners
	bool status = 1;						// status == true => continue to look for points
	while (status == 1){
		start = ntotal - nrecent;
		finish = ntotal;
		nrecent = 0;						// set recent points back to zero for next sweep through
		for (s=start;s<finish;s++){
			// Loop over recent points; look for new points
			x = temp(0,s);
			y = temp(1,s);
			z = temp(2,s);
			// Looop over the directions
			for (p=0;p<26;p++){
				nodx=x+d[p][0];
				if (nodx < 0 ){ nodx = m-1; }		// Periodic BC for x
				if (nodx > m-1 ){ nodx = 0; }
				nody=y+d[p][1];
				if (nody < 0 ){ nody = n-1; }	// Periodic BC for y
				if (nody > n-1 ){ nody = 0; }
				nodz=z+d[p][2];
				if (nodz < 0 ){ nodz = 0; }		// No periodic BC for z
				if (nodz > o-1 ){ nodz = o-1; }
				if ( F(nodx,nody,nodz) > vf && S(nodx,nody,nodz) > vs
					 && indicator(nodx,nody,nodz) == -1 ){
					// Node is a part of the blob - add it to the list
					temp(0,ntotal) = nodx;
					temp(1,ntotal) = nody;
					temp(2,ntotal) = nodz;
					ntotal++;
					nrecent++;
					// Update the indicator map
					indicator(nodx,nody,nodz) = nblobs;
					// Update the min / max for the cube loop
					if ( nodx < imin ){ imin = nodx; }
					if ( nodx > imax ){ imax = nodx; }
					if ( nody < jmin ){ jmin = nody; }
					if ( nody > jmax ){ jmax = nody; }
					if ( nodz < kmin ){ kmin = nodz; }
					if ( nodz > kmax ){ kmax = nodz; }
				}
				else if (F(nodx,nody,nodz) > vf && S(nodx,nody,nodz) > vs
						 && indicator(nodx,nody,nodz) > -1 &&  indicator(nodx,nody,nodz) != nblobs){
					// Some kind of error in algorithm
					printf("Error in blob search algorithm!");
				}
			}

		}
		if ( nrecent == 0){
			status = 0;
		}
	}
	// Use points in temporary storage array to add cubes to the list of blobs
	if ( imin > 0) { imin = imin-1; }
//	if ( imax < m-1) { imax = imax+1; }
	if ( jmin > 0) { jmin = jmin-1; }
//	if ( jmax < n-1) { jmax = jmax+1; }
	if ( kmin > 0) { kmin = kmin-1; }
//	if ( kmax < o-1) { kmax = kmax+1; }
	int i,j,k;
	bool add;
	for (k=kmin;k<kmax;k++){
		for (j=jmin;j<jmax;j++){
			for (i=imin;i<imax;i++){
				// If cube(i,j,k) has any nodes in blob, add it to the list
				// Loop over cube edges
				add = 0;
				for (p=0;p<8;p++){
					nodx = i+cube[p][0];
					nody = j+cube[p][1];
					nodz = k+cube[p][2];
					if ( indicator(nodx,nody,nodz) == nblobs ){
						// Cube corner is in this blob
						add = 1;
					}
				}
				if (add == 1){
					// Add cube to the list
					blobs(0,ncubes) = i;
					blobs(1,ncubes) = j;
					blobs(2,ncubes) = k;
					ncubes++;
					cubes_in_blob++;
					// Loop again to check for overlap
					for (p=0;p<8;p++){
						nodx = i+cube[p][0];
						nody = j+cube[p][1];
						nodz = k+cube[p][2];
						if (indicator(nodx,nody,nodz) > -1 && indicator(nodx,nody,nodz) != nblobs){
							printf("Overlapping cube!");
							cout << i << ", " << j << ", " << k << endl;
						}
					}
				}
			}
		}
	}

	return cubes_in_blob;
}

int main()
{
	int m,n,o;
	int p,i,j,k;

	double vF,vS;
	/* ****************************************************************
			IDENTIFY ALL BLOBS: F > vF, S > vS
	 ****************************************************************** */
	// Find blob domains, number of blobs
	int nblobs = 0;					// number of blobs
	int ncubes = 0;					// total number of nodes in any blob
	int N = (m-1)*(n-1)*(o-1);		// total number of cubes
	IntArray blobs(3,N);	// store indices for blobs (cubes)
	IntArray temp(3,N);		// temporary storage array
	IntArray temp2;			// temporary storage array
	IntArray  b(50);		// number of nodes in each blob

	IntArray indicator(m,n,o);
	DoubleArray F(m,n,o);
	DoubleArray S(m,n,o);

	// Loop over z=0 first -> blobs attached to this end considered "connected" for LB simulation
	i=0;
	int number=0;
	for (k=0;k<1;k++){
		for (j=0;j<n;j++){
			if ( F(i,j,k) > vF ){
				if ( S(i,j,k) > vS ){
					// node i,j,k is in the porespace
					number = number+ComputeBlob(blobs,nblobs,ncubes,indicator,F,S,vF,vS,i,j,k,temp);
				}
			}
	}
	// Specify the blob on the z axis
	b(nblobs) = number;
	nblobs++;
	for (k=0;k<o;k++){
		for (j=0;j<n;j++){
			for (i=1;i<m;i++){
				if ( indicator(i,j,k) == -1 ){
					if ( F(i,j,k) > vF ){
						if ( S(i,j,k) > vS ){
							// node i,j,k is in the porespace
							b(nblobs) = ComputeBlob(blobs,nblobs,ncubes,indicator,F,S,vF,vS,i,j,k,temp);
							nblobs++;
						}
					}
				}
				// Otherwise, this point has already been assigned - ignore

				// Make sure list blob_nodes is large enough
				if ( nblobs > b.length-1){
					b = IncreaseSize(b,b.length);
				}
			}
		}
	}
	// Go over all cubes again -> add any that do not contain nw phase
	bool add=1;			// Set to false if any corners contain nw-phase ( F > vF)
	int cube[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};  // cube corners
	int count_in=0,count_out=0;
	int nodx,nody,nodz;
	for (k=0;k<o-1;k++){
		for (j=0;j<n-1;j++){
			for (i=0;i<m-1;i++){
				// Loop over cube corners
				add=1;				// initialize to true - add unless corner occupied by nw-phase
				for (p=0;p<8;p++){
					nodx=i+cube[p][0];
					nody=j+cube[p][1];
					nodz=k+cube[p][2];
					if ( indicator(nodx,nody,nodz) > -1 ){
						// corner occupied by nw-phase  -> do not add
						add = 0;
					}
				}
				if ( add == 1 ){
					blobs(0,ncubes) = i;
					blobs(1,ncubes) = j;
					blobs(2,ncubes) = k;
					ncubes++;
					count_in++;
				}
				else { count_out++; }
			}
		}
	}
	b(nblobs) = count_in;
	nblobs++;
}

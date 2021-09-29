// Morphological opening routine
#include "common/Array.h"
#include "common/Domain.h"
#include "analysis/runAnalysis.h"

double MorphOpen(DoubleArray &SignDist, signed char *id, std::shared_ptr<Domain> Dm, double VoidFraction, signed char ErodeLabel, signed char ReplaceLabel);
double MorphDrain(DoubleArray &SignDist, signed char *id, std::shared_ptr<Domain> Dm, double VoidFraction);
double MorphGrow(DoubleArray &BoundaryDist, DoubleArray &Dist, Array<char> &id, std::shared_ptr<Domain> Dm, double TargetVol, double WallFactor);

#ifndef MORPHOLOGY_INC
#define MORPHOLOGY_INC
/**
 * \class Morphology
 * @brief 
 * The Morphology class supports morphological operations on complex structures
 * 
 */
class Morphology{
public:    
	/**
	* \brief Create a flow adaptor to operate on the LB model
	* @param Dm       Domain structure 
	*/
	Morphology(std::shared_ptr <Domain> Dm);

	/**
	* \brief Destructor
	*/
	~Morphology();
	
	
	/**
	* \brief Send /recieve function for labels
	*/
	void SendRecv(std::shared_ptr <Domain> Dm, char *id);
	
	/**
	* \brief  Find all sites such that the reach of the signed distance at the site  overlaps with a sub-domain boundary
	*/
	void GetOverlaps(DoubleArray &SignDist);


private:
	IntArray LocalOverlaps;
	IntArray NonLocalOverlaps; 
	
	// List of sites that overlap with processor boundaries
	int *sendOverlap_x, *sendOverlap_y, *sendOverlap_z, *sendOverlap_X, *sendOverlap_Y, *sendOverlap_Z;
	int *sendOverlap_xy, *sendOverlap_yz, *sendOverlap_xz, *sendOverlap_Xy, *sendOverlap_Yz, *sendOverlap_xZ;
	int *sendOverlap_xY, *sendOverlap_yZ, *sendOverlap_Xz, *sendOverlap_XY, *sendOverlap_YZ, *sendOverlap_XZ;
	
	int *recvOverlap_x, *recvOverlap_y, *recvOverlap_z, *recvOverlap_X, *recvOverlap_Y, *recvOverlap_Z;
	int *recvOverlap_xy, *recvOverlap_yz, *recvOverlap_xz, *recvOverlap_Xy, *recvOverlap_Yz, *recvOverlap_xZ;
	int *recvOverlap_xY, *recvOverlap_yZ, *recvOverlap_Xz, *recvOverlap_XY, *recvOverlap_YZ, *recvOverlap_XZ;
	
	// Communication buffers
	signed char *sendID_x, *sendID_y, *sendID_z, *sendID_X, *sendID_Y, *sendID_Z;
	signed char *sendID_xy, *sendID_yz, *sendID_xz, *sendID_Xy, *sendID_Yz, *sendID_xZ;
	signed char *sendID_xY, *sendID_yZ, *sendID_Xz, *sendID_XY, *sendID_YZ, *sendID_XZ;
	signed char *recvID_x, *recvID_y, *recvID_z, *recvID_X, *recvID_Y, *recvID_Z;
	signed char *recvID_xy, *recvID_yz, *recvID_xz, *recvID_Xy, *recvID_Yz, *recvID_xZ;
	signed char *recvID_xY, *recvID_yZ, *recvID_Xz, *recvID_XY, *recvID_YZ, *recvID_XZ;
};

#endif
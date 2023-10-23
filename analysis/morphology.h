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
// Morphological opening routine
#include "common/Array.h"
#include "common/Domain.h"
#include "analysis/runAnalysis.h"

double MorphOpen(DoubleArray &SignDist, signed char *id,
                 std::shared_ptr<Domain> Dm, double VoidFraction,
                 signed char ErodeLabel, signed char ReplaceLabel);
double MorphDrain(DoubleArray &SignDist, signed char *id,
                  std::shared_ptr<Domain> Dm, double VoidFraction,
                  double InitialRadius);
double MorphGrow(DoubleArray &BoundaryDist, DoubleArray &Dist, Array<char> &id,
                 std::shared_ptr<Domain> Dm, double TargetVol,
                 double WallFactor);

#ifndef MORPHOLOGY_INC
#define MORPHOLOGY_INC
/**
 * \class Morphology
 * @brief 
 * The Morphology class supports morphological operations on complex structures
 * 
 */
class Morphology {
public:
    /**
	* \brief Create a flow adaptor to operate on the LB model
	*/
    Morphology();

    /**
	* \brief Destructor
	*/
    ~Morphology();

    /**
	* \brief Initialize morphology structure from distance map
	* @param Dm       Domain structure 
	* @param Distance   Signed distance to boundary of structure
	*/
    void Initialize(std::shared_ptr<Domain> Dm, DoubleArray &Distance);

    /**
	* \brief  Find all sites such that the reach of the signed distance at the site  overlaps with a sub-domain boundary
	* @param Dm    Domain structure 
	* @param id    image labels
	* @param ErodeLabel label to erode based on morphological operation
	* @param NewLabel   label to assign based on morphological operation
	*/
    int GetOverlaps(std::shared_ptr<Domain> Dm, signed char *id,
                    const signed char ErodeLabel, const signed char NewLabel);

    /*
     *  data structures to store non-local morphological information 
    */
    std::vector<int> xShift, yShift, zShift;
    std::vector<int> sendID;
    std::vector<double> morphRadius;
    std::vector<unsigned char> localID;
    std::vector<unsigned char> nonlocalID;

private:
    int sendtag, recvtag;

    //......................................................................................
    int sendCount, recvCount;
    //......................................................................................
    int sendOffset_x, sendOffset_y, sendOffset_z, sendOffset_X, sendOffset_Y,
        sendOffset_Z;
    int sendOffset_xy, sendOffset_yz, sendOffset_xz, sendOffset_Xy,
        sendOffset_Yz, sendOffset_xZ;
    int sendOffset_xY, sendOffset_yZ, sendOffset_Xz, sendOffset_XY,
        sendOffset_YZ, sendOffset_XZ;
    int sendOffset_xyz, sendOffset_XYZ, sendOffset_xYz, sendOffset_XyZ;
    int sendOffset_Xyz, sendOffset_xYZ, sendOffset_xyZ, sendOffset_XYz;
    //......................................................................................
    int recvOffset_x, recvOffset_y, recvOffset_z, recvOffset_X, recvOffset_Y,
        recvOffset_Z;
    int recvOffset_xy, recvOffset_yz, recvOffset_xz, recvOffset_Xy,
        recvOffset_Yz, recvOffset_xZ;
    int recvOffset_xY, recvOffset_yZ, recvOffset_Xz, recvOffset_XY,
        recvOffset_YZ, recvOffset_XZ;
    int recvOffset_xyz, recvOffset_XYZ, recvOffset_xYz, recvOffset_XyZ;
    int recvOffset_Xyz, recvOffset_xYZ, recvOffset_xyZ, recvOffset_XYz;
    //......................................................................................
    int sendCount_x, sendCount_y, sendCount_z, sendCount_X, sendCount_Y,
        sendCount_Z;
    int sendCount_xy, sendCount_yz, sendCount_xz, sendCount_Xy, sendCount_Yz,
        sendCount_xZ;
    int sendCount_xY, sendCount_yZ, sendCount_Xz, sendCount_XY, sendCount_YZ,
        sendCount_XZ;
    int sendCount_xyz, sendCount_XYZ, sendCount_xYz, sendCount_XyZ;
    int sendCount_Xyz, sendCount_xYZ, sendCount_xyZ, sendCount_XYz;
    //......................................................................................
    int recvCount_x, recvCount_y, recvCount_z, recvCount_X, recvCount_Y,
        recvCount_Z;
    int recvCount_xy, recvCount_yz, recvCount_xz, recvCount_Xy, recvCount_Yz,
        recvCount_xZ;
    int recvCount_xY, recvCount_yZ, recvCount_Xz, recvCount_XY, recvCount_YZ,
        recvCount_XZ;
    int recvCount_xyz, recvCount_XYZ, recvCount_xYz, recvCount_XyZ;
    int recvCount_Xyz, recvCount_xYZ, recvCount_xyZ, recvCount_XYz;
    //......................................................................................
    std::vector<char> sendList;
    std::vector<char> recvList;
    //......................................................................................

    // Communication buffers
    signed char *sendID_x, *sendID_y, *sendID_z, *sendID_X, *sendID_Y,
        *sendID_Z;
    signed char *sendID_xy, *sendID_yz, *sendID_xz, *sendID_Xy, *sendID_Yz,
        *sendID_xZ;
    signed char *sendID_xY, *sendID_yZ, *sendID_Xz, *sendID_XY, *sendID_YZ,
        *sendID_XZ;
    signed char *recvID_x, *recvID_y, *recvID_z, *recvID_X, *recvID_Y,
        *recvID_Z;
    signed char *recvID_xy, *recvID_yz, *recvID_xz, *recvID_Xy, *recvID_Yz,
        *recvID_xZ;
    signed char *recvID_xY, *recvID_yZ, *recvID_Xz, *recvID_XY, *recvID_YZ,
        *recvID_XZ;
};

#endif
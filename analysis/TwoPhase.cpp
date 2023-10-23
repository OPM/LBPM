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
#include "analysis/TwoPhase.h"

#include "analysis/pmmc.h"
#include "analysis/analysis.h"
#include "common/Domain.h"
#include "common/Communication.h"
#include "common/Utilities.h"
#include "common/MPI.h"
#include "IO/MeshDatabase.h"
#include "IO/Reader.h"
#include "IO/Writer.h"
#include "analysis/filters.h"

#include <memory>

#define BLOB_AVG_COUNT 35

// Array access for averages defined by the following
#define VOL 0
#define TRIMVOL 1
#define PRS 2
#define AWN 3
#define AWS 4
#define ANS 5
#define LWNS 6
#define JWN 7
#define KWN 8
#define CWNS 9
#define KNWNS 10
#define KGWNS 11
#define VX 12
#define VY 13
#define VZ 14
#define VSQ 15
#define VWNX 16
#define VWNY 17
#define VWNZ 18
#define VWNSX 19
#define VWNSY 20
#define VWNSZ 21
#define GWNXX 22
#define GWNYY 23
#define GWNZZ 24
#define GWNXY 25
#define GWNXZ 26
#define GWNYZ 27
#define TRAWN 28
#define TRJWN 29
#define CMX 30
#define CMY 31
#define CMZ 32
#define EULER 33
#define INTCURV 34

#define PI 3.14159265359

// Constructor
TwoPhase::TwoPhase(std::shared_ptr<Domain> dm)
    : n_nw_pts(0), n_ns_pts(0), n_ws_pts(0), n_nws_pts(0), n_local_sol_pts(0),
      n_local_nws_pts(0), n_nw_tris(0), n_ns_tris(0), n_ws_tris(0),
      n_nws_seg(0), n_local_sol_tris(0), nc(0), kstart(0), kfinish(0),
      fluid_isovalue(0), solid_isovalue(0), Volume(0), TIMELOG(NULL),
      NWPLOG(NULL), WPLOG(NULL), Dm(dm), NumberComponents_WP(0),
      NumberComponents_NWP(0), trimdist(0), porosity(0), poreVol(0), awn(0),
      ans(0), aws(0), lwns(0), wp_volume(0), nwp_volume(0), As(0), dummy(0),
      vol_w(0), vol_n(0), sat_w(0), sat_w_previous(0), pan(0), paw(0),
      pan_global(0), paw_global(0), vol_w_global(0), vol_n_global(0),
      awn_global(0), ans_global(0), aws_global(0), lwns_global(0), efawns(0),
      efawns_global(0), Jwn(0), Jwn_global(0), Kwn(0), Kwn_global(0), KNwns(0),
      KNwns_global(0), KGwns(0), KGwns_global(0), trawn(0), trawn_global(0),
      trJwn(0), trJwn_global(0), trRwn(0), trRwn_global(0),
      nwp_volume_global(0), wp_volume_global(0), As_global(0), wwndnw_global(0),
      wwnsdnwn_global(0), Jwnwwndnw_global(0), dEs(0), dAwn(0), dAns(0) {
    Nx = dm->Nx;
    Ny = dm->Ny;
    Nz = dm->Nz;
    Volume = (Nx - 2) * (Ny - 2) * (Nz - 2) * Dm->nprocx() * Dm->nprocy() *
             Dm->nprocz() * 1.0;

    TempID = new char[Nx * Ny * Nz];

    wet_morph = std::shared_ptr<Minkowski>(new Minkowski(Dm));
    nonwet_morph = std::shared_ptr<Minkowski>(new Minkowski(Dm));

    // Global arrays
    PhaseID.resize(Nx, Ny, Nz);
    PhaseID.fill(0);
    Label_WP.resize(Nx, Ny, Nz);
    Label_WP.fill(0);
    Label_NWP.resize(Nx, Ny, Nz);
    Label_NWP.fill(0);
    SDn.resize(Nx, Ny, Nz);
    SDn.fill(0);
    SDs.resize(Nx, Ny, Nz);
    SDs.fill(0);
    Phase.resize(Nx, Ny, Nz);
    Phase.fill(0);
    Press.resize(Nx, Ny, Nz);
    Press.fill(0);
    dPdt.resize(Nx, Ny, Nz);
    dPdt.fill(0);
    MeanCurvature.resize(Nx, Ny, Nz);
    MeanCurvature.fill(0);
    GaussCurvature.resize(Nx, Ny, Nz);
    GaussCurvature.fill(0);
    SDs_x.resize(Nx, Ny, Nz);
    SDs_x.fill(0); // Gradient of the signed distance
    SDs_y.resize(Nx, Ny, Nz);
    SDs_y.fill(0);
    SDs_z.resize(Nx, Ny, Nz);
    SDs_z.fill(0);
    SDn_x.resize(Nx, Ny, Nz);
    SDn_x.fill(0); // Gradient of the signed distance
    SDn_y.resize(Nx, Ny, Nz);
    SDn_y.fill(0);
    SDn_z.resize(Nx, Ny, Nz);
    SDn_z.fill(0);
    DelPhi.resize(Nx, Ny, Nz);
    DelPhi.fill(0);
    Phase_tplus.resize(Nx, Ny, Nz);
    Phase_tplus.fill(0);
    Phase_tminus.resize(Nx, Ny, Nz);
    Phase_tminus.fill(0);
    Vel_x.resize(Nx, Ny, Nz);
    Vel_x.fill(0); // Gradient of the phase indicator field
    Vel_y.resize(Nx, Ny, Nz);
    Vel_y.fill(0);
    Vel_z.resize(Nx, Ny, Nz);
    Vel_z.fill(0);
    //.........................................
    // Allocate cube storage space
    CubeValues.resize(2, 2, 2);
    nw_tris.resize(3, 20);
    ns_tris.resize(3, 20);
    ws_tris.resize(3, 20);
    nws_seg.resize(2, 20);
    local_sol_tris.resize(3, 18);
    nw_pts = DTMutableList<Point>(20);
    ns_pts = DTMutableList<Point>(20);
    ws_pts = DTMutableList<Point>(20);
    nws_pts = DTMutableList<Point>(20);
    local_nws_pts = DTMutableList<Point>(20);
    local_sol_pts = DTMutableList<Point>(20);
    tmp = DTMutableList<Point>(20);
    //.........................................
    Values.resize(20);
    DistanceValues.resize(20);
    KGwns_values.resize(20);
    KNwns_values.resize(20);
    InterfaceSpeed.resize(20);
    NormalVector.resize(60);
    //.........................................
    van.resize(3);
    vaw.resize(3);
    vawn.resize(3);
    vawns.resize(3);
    Gwn.resize(6);
    Gns.resize(6);
    Gws.resize(6);
    van_global.resize(3);
    vaw_global.resize(3);
    vawn_global.resize(3);
    vawns_global.resize(3);
    Gwn_global.resize(6);
    Gns_global.resize(6);
    Gws_global.resize(6);
    //.........................................
    if (Dm->rank() == 0) {
        TIMELOG = fopen("timelog.tcat", "a+");
        if (fseek(TIMELOG, 0, SEEK_SET) == fseek(TIMELOG, 0, SEEK_CUR)) {
            // If timelog is empty, write a short header to list the averages
            //fprintf(TIMELOG,"--------------------------------------------------------------------------------------\n");
            fprintf(
                TIMELOG,
                "time rn rw nun nuw Fx Fy Fz iftwn "); // Timestep, Change in Surface Energy
            fprintf(TIMELOG, "sw pw pn awn ans aws Jwn Kwn lwns cwns KNwns "
                             "KGwns "); // Scalar averages
            fprintf(TIMELOG,
                    "vawx vawy vawz vanx vany vanz "); // Velocity averages
            fprintf(TIMELOG, "vawnx vawny vawnz vawnsx vawnsy vawnsz ");
            fprintf(
                TIMELOG,
                "Gwnxx Gwnyy Gwnzz Gwnxy Gwnxz Gwnyz "); // Orientation tensors
            fprintf(TIMELOG, "Gwsxx Gwsyy Gwszz Gwsxy Gwsxz Gwsyz ");
            fprintf(TIMELOG, "Gnsxx Gnsyy Gnszz Gnsxy Gnsxz Gnsyz ");
            fprintf(TIMELOG, "trawn trJwn trRwn "); //trimmed curvature,
            fprintf(TIMELOG,
                    "wwndnw wwnsdnwn Jwnwwndnw "); //kinematic quantities,
            fprintf(TIMELOG, "Vw Aw Jw Xw ");      //miknowski measures,
            fprintf(TIMELOG, "Vn An Jn Xn\n");     //miknowski measures,
            //			fprintf(TIMELOG,"Euler Kn Jn An\n"); 			//miknowski measures,
        }

        NWPLOG = fopen("components.NWP.tcat", "a+");
        fprintf(NWPLOG, "time label vol pn awn ans Jwn Kwn lwns cwns ");
        fprintf(NWPLOG, "vx vy vz vwnx vwny vwnz vwnsx vwnsy vwnsz vsq ");
        fprintf(NWPLOG, "Gwnxx Gwnyy Gwnzz Gwnxy Gwnxz Gwnyz Cx Cy Cz trawn "
                        "trJwn Kn Euler\n");

        WPLOG = fopen("components.WP.tcat", "a+");
        fprintf(WPLOG, "time label vol pw awn ans Jwn Kwn lwns cwns ");
        fprintf(WPLOG, "vx vy vz vwnx vwny vwnz vwnsx vwnsy vwnsz vsq ");
        fprintf(WPLOG,
                "Gwnxx Gwnyy Gwnzz Gwnxy Gwnxz Gwnyz Cx Cy Cz trawn trJwn\n");
    } else {
        char LocalRankString[8];
        sprintf(LocalRankString, "%05d", Dm->rank());
        char LocalRankFilename[40];
        sprintf(LocalRankFilename, "%s%s", "timelog.tcat.", LocalRankString);
        TIMELOG = fopen(LocalRankFilename, "a+");
        //fprintf(TIMELOG,"--------------------------------------------------------------------------------------\n");
        fprintf(TIMELOG, "time rn rw nun nuw Fx Fy Fz iftwn ");
        ; // Timestep,
        fprintf(
            TIMELOG,
            "sw pw pn awn ans aws Jwn Kwn lwns cwns KNwns KGwns "); // Scalar averages
        fprintf(TIMELOG, "vawx vawy vawz vanx vany vanz "); // Velocity averages
        fprintf(TIMELOG, "vawnx vawny vawnz vawnsx vawnsy vawnsz ");
        fprintf(TIMELOG,
                "Gwnxx Gwnyy Gwnzz Gwnxy Gwnxz Gwnyz "); // Orientation tensors
        fprintf(TIMELOG, "Gwsxx Gwsyy Gwszz Gwsxy Gwsxz Gwsyz ");
        fprintf(TIMELOG, "Gnsxx Gnsyy Gnszz Gnsxy Gnsxz Gnsyz ");
        fprintf(TIMELOG, "trawn trJwn trRwn ");         //trimmed curvature,
        fprintf(TIMELOG, "wwndnw wwnsdnwn Jwnwwndnw "); //kinematic quantities,
        fprintf(TIMELOG, "Vw Aw Jw Xw ");               //miknowski measures,
        fprintf(TIMELOG, "Vn An Jn Xn\n");              //miknowski measures,
        //	fprintf(TIMELOG,"Euler Kn Jn An\n"); 			//miknowski measures,
    }
}

// Destructor
TwoPhase::~TwoPhase() {
    delete[] TempID;
    if (TIMELOG != NULL) {
        fclose(TIMELOG);
    }
    if (NWPLOG != NULL) {
        fclose(NWPLOG);
    }
    if (WPLOG != NULL) {
        fclose(WPLOG);
    }
}

void TwoPhase::ColorToSignedDistance(double Beta, DoubleArray &ColorData,
                                     DoubleArray &DistData) {
    NULL_USE(Beta);
    /*double factor,temp,value;
		factor=0.5/Beta;
	// Initialize to -1,1 (segmentation)
	for (int k=0; k<Nz; k++){
		for (int j=0; j<Ny; j++){
			for (int i=0; i<Nx; i++){
				value = ColorData(i,j,k);

				// Set phase ID field for non-wetting phase
				// here NWP is tagged with 1
				// all other phases tagged with 0
				int n = k*Nx*Ny+j*Nx+i;
				if (value > 0)	TempID[n] = 1;
				else		    TempID[n] = 0;

				// Distance threshhold 
				// temp -- distance based on analytical form McClure, Prins et al, Comp. Phys. Comm.
				//  distance should be negative outside the NWP
				//  distance should be positive inside of the NWP
				temp = factor*log((1.0+value)/(1.0-value));
				if (value > 0.8) DistData(i,j,k) = 2.94*factor;
				else if (value < -0.8) DistData(i,j,k) = -2.94*factor;
				else DistData(i,j,k) = temp;

				// Basic threshold
				// distance to the NWP
				// negative inside NWP, positive outside
				//if (value > 0) DistData(i,j,k) = -0.5;
				//else DistData(i,j,k) = 0.5;

				// Initialize directly
				// DistData(i,j,k) = value;
			}
		}
	}

	CalcDist(DistData,TempID,Dm);

	for (int k=0; k<Nz; k++){
		for (int j=0; j<Ny; j++){
			for (int i=0; i<Nx; i++){
				DistData(i,j,k) += 1.0;
			}
		}
	}	
	*/
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                DistData(i, j, k) = ColorData(i, j, k);
            }
        }
    }
}

void TwoPhase::ComputeDelPhi() {
    int i, j, k;
    double fx, fy, fz;

    Dm->CommunicateMeshHalo(Phase);
    for (k = 1; k < Nz - 1; k++) {
        for (j = 1; j < Ny - 1; j++) {
            for (i = 1; i < Nx - 1; i++) {
                // Compute all of the derivatives using finite differences
                fx = 0.5 * (Phase(i + 1, j, k) - Phase(i - 1, j, k));
                fy = 0.5 * (Phase(i, j + 1, k) - Phase(i, j - 1, k));
                fz = 0.5 * (Phase(i, j, k + 1) - Phase(i, j, k - 1));
                DelPhi(i, j, k) = sqrt(fx * fx + fy * fy + fz * fz);
            }
        }
    }
}

void TwoPhase::Initialize() {
    trimdist = 1.0;
    fluid_isovalue = solid_isovalue = 0.0;
    // Initialize the averaged quantities
    awn = aws = ans = lwns = 0.0;
    nwp_volume = wp_volume = 0.0;
    As = 0.0;
    pan = paw = 0.0;
    vaw(0) = vaw(1) = vaw(2) = 0.0;
    van(0) = van(1) = van(2) = 0.0;
    vawn(0) = vawn(1) = vawn(2) = 0.0;
    vawns(0) = vawns(1) = vawns(2) = 0.0;
    Gwn(0) = Gwn(1) = Gwn(2) = 0.0;
    Gwn(3) = Gwn(4) = Gwn(5) = 0.0;
    Gws(0) = Gws(1) = Gws(2) = 0.0;
    Gws(3) = Gws(4) = Gws(5) = 0.0;
    Gns(0) = Gns(1) = Gns(2) = 0.0;
    Gns(3) = Gns(4) = Gns(5) = 0.0;
    vol_w = vol_n = 0.0;
    KGwns = KNwns = 0.0;
    Jwn = Kwn = efawns = 0.0;
    trJwn = trawn = trRwn = 0.0;
    euler = Jn = An = Kn = 0.0;
    wwndnw = 0.0;
    wwnsdnwn = 0.0;
    Jwnwwndnw = 0.0;
    Xwn = Xns = Xws = 0.0;
}

void TwoPhase::SetParams(double rhoA, double rhoB, double tauA, double tauB,
                         double force_x, double force_y, double force_z,
                         double alpha) {
    Fx = force_x;
    Fy = force_y;
    Fz = force_z;
    rho_n = rhoA;
    rho_w = rhoB;
    nu_n = (tauA - 0.5) / 3.f;
    nu_w = (tauB - 0.5) / 3.f;
    gamma_wn = 5.796 * alpha;
}

/*
void TwoPhase::SetupCubes(Domain &Dm)
{
	int i,j,k;
	kstart = 1;
	kfinish = Nz-1;
	if (Dm->BoundaryCondition !=0 && Dm->kproc==0)			kstart = 4;
	if (Dm->BoundaryCondition !=0 && Dm->kproc==Dm->nprocz-1)	kfinish = Nz-4;
	nc=0;
	for (k=kstart; k<kfinish; k++){
		for (j=1; j<Ny-1; j++){
			for (i=1; i<Nx-1; i++){
				cubeList(0,nc) = i;
				cubeList(1,nc) = j;
				cubeList(2,nc) = k;
				nc++;
			}
		}
	}
	ncubes = nc;
}
*/

void TwoPhase::UpdateSolid() {
    Dm->CommunicateMeshHalo(SDs);
    //...........................................................................
    // Gradient of the Signed Distance function
    //...........................................................................
    pmmc_MeshGradient(SDs, SDs_x, SDs_y, SDs_z, Nx, Ny, Nz);
    //...........................................................................
    Dm->CommunicateMeshHalo(SDs_x);
    //...........................................................................
    Dm->CommunicateMeshHalo(SDs_y);
    //...........................................................................
    Dm->CommunicateMeshHalo(SDs_z);
    //...........................................................................
}

void TwoPhase::UpdateMeshValues() {
    int i, j, k, n;
    fillHalo<double> fillData(Dm->Comm, Dm->rank_info, {Nx - 2, Ny - 2, Nz - 2},
                              {1, 1, 1}, 0, 1);

    //...........................................................................
    //Dm->CommunicateMeshHalo(SDn);
    fillData.fill(SDn);
    //...........................................................................
    // Compute the gradients of the phase indicator and signed distance fields
    pmmc_MeshGradient(SDn, SDn_x, SDn_y, SDn_z, Nx, Ny, Nz);
    //...........................................................................
    // Gradient of the phase indicator field
    fillData.fill(SDn_x);
    fillData.fill(SDn_y);
    fillData.fill(SDn_z);
    fillData.fill(SDs);
    //...........................................................................
    //Dm->CommunicateMeshHalo(SDn_x);
    //...........................................................................
    //Dm->CommunicateMeshHalo(SDn_y);
    //...........................................................................
    //Dm->CommunicateMeshHalo(SDn_z);
    //...........................................................................
    //Dm->CommunicateMeshHalo(SDs);
    pmmc_MeshGradient(SDs, SDs_x, SDs_y, SDs_z, Nx, Ny, Nz);
    //...........................................................................
    fillData.fill(SDs_x);
    fillData.fill(SDs_y);
    fillData.fill(SDs_z);
    //Dm->CommunicateMeshHalo(SDs_x);
    //...........................................................................
    //Dm->CommunicateMeshHalo(SDs_y);
    //...........................................................................
    //Dm->CommunicateMeshHalo(SDs_z);
    //...........................................................................
    // Compute the mesh curvature of the phase indicator field
    pmmc_MeshCurvature(SDn, MeanCurvature, GaussCurvature, Nx, Ny, Nz);
    //...........................................................................
    // Update the time derivative of non-dimensional density field
    // Map Phase_tplus and Phase_tminus
    for (int n = 0; n < Nx * Ny * Nz; n++)
        dPdt(n) = 0.125 * (Phase_tplus(n) - Phase_tminus(n));
    //...........................................................................
    Dm->CommunicateMeshHalo(Press);
    //...........................................................................
    Dm->CommunicateMeshHalo(Vel_x);
    //...........................................................................
    Dm->CommunicateMeshHalo(Vel_y);
    //...........................................................................
    Dm->CommunicateMeshHalo(Vel_z);
    //...........................................................................
    Dm->CommunicateMeshHalo(MeanCurvature);
    //...........................................................................
    Dm->CommunicateMeshHalo(GaussCurvature);
    //...........................................................................
    Dm->CommunicateMeshHalo(DelPhi);
    //...........................................................................
    // Initializing the blob ID
    for (k = 0; k < Nz; k++) {
        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                n = k * Nx * Ny + j * Nx + i;
                if (!(Dm->id[n] > 0)) {
                    // Solid phase
                    PhaseID(i, j, k) = 0;
                } else if (SDn(i, j, k) < 0.0) {
                    // wetting phase
                    PhaseID(i, j, k) = 2;
                } else {
                    // non-wetting phase
                    PhaseID(i, j, k) = 1;
                }
            }
        }
    }
}
void TwoPhase::ComputeLocal() {
    int i, j, k, n, imin, jmin, kmin, kmax;
    int cube[8][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0},
                      {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};

    // If external boundary conditions are set, do not average over the inlet
    kmin = 1;
    kmax = Nz - 1;
    if (Dm->BoundaryCondition > 0 && Dm->kproc() == 0)
        kmin = 4;
    if (Dm->BoundaryCondition > 0 && Dm->kproc() == Dm->nprocz() - 1)
        kmax = Nz - 4;

    imin = jmin = 1;
    // If inlet layers exist use these as default
    if (Dm->inlet_layers_x > 0)
        imin = Dm->inlet_layers_x;
    if (Dm->inlet_layers_y > 0)
        jmin = Dm->inlet_layers_y;
    if (Dm->inlet_layers_z > 0)
        kmin = Dm->inlet_layers_z;

    for (k = kmin; k < kmax; k++) {
        for (j = jmin; j < Ny - 1; j++) {
            for (i = imin; i < Nx - 1; i++) {
                //...........................................................................
                n_nw_pts = n_ns_pts = n_ws_pts = n_nws_pts = n_local_sol_pts =
                    n_local_nws_pts = 0;
                n_nw_tris = n_ns_tris = n_ws_tris = n_nws_seg =
                    n_local_sol_tris = 0;
                //...........................................................................
                // Compute volume averages
                for (int p = 0; p < 8; p++) {
                    n = i + cube[p][0] + (j + cube[p][1]) * Nx +
                        (k + cube[p][2]) * Nx * Ny;
                    if (Dm->id[n] > 0) {
                        // 1-D index for this cube corner
                        // compute the norm of the gradient of the phase indicator field
                        // Compute the non-wetting phase volume contribution
                        if (Phase(i + cube[p][0], j + cube[p][1],
                                  k + cube[p][2]) > 0) {
                            nwp_volume += 0.125;
                            // velocity
                            van(0) += 0.125 * Vel_x(n);
                            van(1) += 0.125 * Vel_y(n);
                            van(2) += 0.125 * Vel_z(n);
                            // volume the excludes the interfacial region
                            if (DelPhi(n) < 1e-4) {
                                vol_n += 0.125;
                                // pressure
                                pan += 0.125 * Press(n);
                            }
                        } else {
                            wp_volume += 0.125;
                            // velocity
                            vaw(0) += 0.125 * Vel_x(n);
                            vaw(1) += 0.125 * Vel_y(n);
                            vaw(2) += 0.125 * Vel_z(n);
                            if (DelPhi(n) < 1e-4) {
                                // volume the excludes the interfacial region
                                vol_w += 0.125;
                                // pressure
                                paw += 0.125 * Press(n);
                            }
                        }
                    }
                }

                //...........................................................................
                // Construct the interfaces and common curve
                pmmc_ConstructLocalCube(
                    SDs, SDn, solid_isovalue, fluid_isovalue, nw_pts, nw_tris,
                    Values, ns_pts, ns_tris, ws_pts, ws_tris, local_nws_pts,
                    nws_pts, nws_seg, local_sol_pts, local_sol_tris,
                    n_local_sol_tris, n_local_sol_pts, n_nw_pts, n_nw_tris,
                    n_ws_pts, n_ws_tris, n_ns_tris, n_ns_pts, n_local_nws_pts,
                    n_nws_pts, n_nws_seg, i, j, k, Nx, Ny, Nz);

                // wn interface averages
                if (n_nw_pts > 0) {
                    awn += pmmc_CubeSurfaceOrientation(Gwn, nw_pts, nw_tris,
                                                       n_nw_tris);
                    Kwn += pmmc_CubeSurfaceInterpValue(
                        CubeValues, GaussCurvature, nw_pts, nw_tris, Values, i,
                        j, k, n_nw_pts, n_nw_tris);

                    Jwn += pmmc_CubeSurfaceInterpValue(
                        CubeValues, MeanCurvature, nw_pts, nw_tris, Values, i,
                        j, k, n_nw_pts, n_nw_tris);

                    // Compute the normal speed of the interface
                    wwndnw += pmmc_InterfaceSpeed(
                        dPdt, SDn_x, SDn_y, SDn_z, CubeValues, nw_pts, nw_tris,
                        NormalVector, InterfaceSpeed, vawn, i, j, k, n_nw_pts,
                        n_nw_tris);

                    //for (int p=0; p <n_nw_tris; p++) wwndnw += InterfaceSpeed(p);
                    for (int p = 0; p < n_nw_tris; p++)
                        Jwnwwndnw += InterfaceSpeed(p) * Values(p);

                    // Integrate the trimmed mean curvature (hard-coded to use a distance of 4 pixels)
                    pmmc_CubeTrimSurfaceInterpValues(
                        CubeValues, MeanCurvature, SDs, nw_pts, nw_tris, Values,
                        DistanceValues, i, j, k, n_nw_pts, n_nw_tris, trimdist,
                        dummy, trJwn);

                    pmmc_CubeTrimSurfaceInterpInverseValues(
                        CubeValues, MeanCurvature, SDs, nw_pts, nw_tris, Values,
                        DistanceValues, i, j, k, n_nw_pts, n_nw_tris, trimdist,
                        dummy, trRwn);
                }
                // wns common curve averages
                if (n_local_nws_pts > 0) {
                    efawns += pmmc_CubeContactAngle(
                        CubeValues, Values, SDn_x, SDn_y, SDn_z, SDs_x, SDs_y,
                        SDs_z, local_nws_pts, i, j, k, n_local_nws_pts);

                    wwnsdnwn += pmmc_CommonCurveSpeed(
                        CubeValues, dPdt, vawns, SDn_x, SDn_y, SDn_z, SDs_x,
                        SDs_y, SDs_z, local_nws_pts, i, j, k, n_local_nws_pts);

                    pmmc_CurveCurvature(SDn, SDs, SDn_x, SDn_y, SDn_z, SDs_x,
                                        SDs_y, SDs_z, KNwns_values,
                                        KGwns_values, KNwns, KGwns, nws_pts,
                                        n_nws_pts, i, j, k);

                    lwns +=
                        pmmc_CubeCurveLength(local_nws_pts, n_local_nws_pts);
                }

                // Solid interface averagees
                if (n_local_sol_tris > 0) {
                    As += pmmc_CubeSurfaceArea(local_sol_pts, local_sol_tris,
                                               n_local_sol_tris);

                    // Compute the surface orientation and the interfacial area
                    ans += pmmc_CubeSurfaceOrientation(Gns, ns_pts, ns_tris,
                                                       n_ns_tris);
                    aws += pmmc_CubeSurfaceOrientation(Gws, ws_pts, ws_tris,
                                                       n_ws_tris);
                }
                //...........................................................................
                // Compute the integral curvature of the non-wetting phase

                n_nw_pts = n_nw_tris = 0;
                // Compute the non-wetting phase surface and associated area
                An +=
                    geomavg_MarchingCubes(SDn, fluid_isovalue, i, j, k, nw_pts,
                                          n_nw_pts, nw_tris, n_nw_tris);
                // Compute the integral of mean curvature
                if (n_nw_pts > 0) {
                    pmmc_CubeTrimSurfaceInterpValues(
                        CubeValues, MeanCurvature, SDs, nw_pts, nw_tris, Values,
                        DistanceValues, i, j, k, n_nw_pts, n_nw_tris, trimdist,
                        trawn, dummy);
                }

                Jn += pmmc_CubeSurfaceInterpValue(CubeValues, MeanCurvature,
                                                  nw_pts, nw_tris, Values, i, j,
                                                  k, n_nw_pts, n_nw_tris);
                // Compute Euler characteristic from integral of gaussian curvature
                Kn += pmmc_CubeSurfaceInterpValue(CubeValues, GaussCurvature,
                                                  nw_pts, nw_tris, Values, i, j,
                                                  k, n_nw_pts, n_nw_tris);

                euler += geomavg_EulerCharacteristic(nw_pts, nw_tris, n_nw_pts,
                                                     n_nw_tris, i, j, k);
            }
        }
    }

    Array<char> phase_label(Nx, Ny, Nz);
    Array<double> phase_distance(Nx, Ny, Nz);
    // Analyze the wetting fluid
    for (k = 0; k < Nz; k++) {
        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                n = k * Nx * Ny + j * Nx + i;
                if (!(Dm->id[n] > 0)) {
                    // Solid phase
                    phase_label(i, j, k) = 1;
                } else if (SDn(i, j, k) < 0.0) {
                    // wetting phase
                    phase_label(i, j, k) = 0;
                } else {
                    // non-wetting phase
                    phase_label(i, j, k) = 1;
                }
                phase_distance(i, j, k) =
                    2.0 * double(phase_label(i, j, k)) - 1.0;
            }
        }
    }
    CalcDist(phase_distance, phase_label, *Dm);
    wet_morph->ComputeScalar(phase_distance, 0.f);
    //printf("generating distance at rank=%i \n",Dm->rank());
    // Analyze the wetting fluid
    for (k = 0; k < Nz; k++) {
        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                n = k * Nx * Ny + j * Nx + i;
                if (!(Dm->id[n] > 0)) {
                    // Solid phase
                    phase_label(i, j, k) = 1;
                } else if (SDn(i, j, k) < 0.0) {
                    // wetting phase
                    phase_label(i, j, k) = 1;
                } else {
                    // non-wetting phase
                    phase_label(i, j, k) = 0;
                }
                phase_distance(i, j, k) =
                    2.0 * double(phase_label(i, j, k)) - 1.0;
            }
        }
    }
    //printf("calculate distance at rank=%i \n",Dm->rank());
    CalcDist(phase_distance, phase_label, *Dm);
    //printf("morphological analysis at rank=%i \n",Dm->rank());
    nonwet_morph->ComputeScalar(phase_distance, 0.f);
    //printf("rank=%i completed \n",Dm->rank());
}
void TwoPhase::ComputeStatic() {
    int i, j, k, n, imin, jmin, kmin, kmax;
    int cube[8][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0},
                      {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};

    kmin = 1;
    kmax = Nz - 1;
    imin = jmin = 1;

    /* set fluid isovalue to "grow" NWP for contact angle measurement */
    fluid_isovalue = -1.0;

    string FILENAME = "ContactAngle";

    char LocalRankString[8];
    char LocalRankFilename[40];
    sprintf(LocalRankString, "%05d", Dm->rank());
    sprintf(LocalRankFilename, "%s%s%s", "ContactAngle.", LocalRankString,
            ".csv");
    FILE *ANGLES = fopen(LocalRankFilename, "a+");
    fprintf(ANGLES, "x y z angle\n");

    for (k = kmin; k < kmax; k++) {
        for (j = jmin; j < Ny - 1; j++) {
            for (i = imin; i < Nx - 1; i++) {
                //...........................................................................
                n_nw_pts = n_ns_pts = n_ws_pts = n_nws_pts = n_local_sol_pts =
                    n_local_nws_pts = 0;
                n_nw_tris = n_ns_tris = n_ws_tris = n_nws_seg =
                    n_local_sol_tris = 0;
                //...........................................................................
                // Compute volume averages
                for (int p = 0; p < 8; p++) {
                    n = i + cube[p][0] + (j + cube[p][1]) * Nx +
                        (k + cube[p][2]) * Nx * Ny;
                    if (Dm->id[n] > 0) {
                        // 1-D index for this cube corner
                        // compute the norm of the gradient of the phase indicator field
                        // Compute the non-wetting phase volume contribution
                        if (Phase(i + cube[p][0], j + cube[p][1],
                                  k + cube[p][2]) > 0) {
                            nwp_volume += 0.125;
                        } else {
                            wp_volume += 0.125;
                        }
                    }
                }

                //...........................................................................
                // Construct the interfaces and common curve
                pmmc_ConstructLocalCube(
                    SDs, SDn, solid_isovalue, fluid_isovalue, nw_pts, nw_tris,
                    Values, ns_pts, ns_tris, ws_pts, ws_tris, local_nws_pts,
                    nws_pts, nws_seg, local_sol_pts, local_sol_tris,
                    n_local_sol_tris, n_local_sol_pts, n_nw_pts, n_nw_tris,
                    n_ws_pts, n_ws_tris, n_ns_tris, n_ns_pts, n_local_nws_pts,
                    n_nws_pts, n_nws_seg, i, j, k, Nx, Ny, Nz);

                // wn interface averages
                if (n_nw_pts > 0) {
                    awn += pmmc_CubeSurfaceOrientation(Gwn, nw_pts, nw_tris,
                                                       n_nw_tris);
                    Kwn += pmmc_CubeSurfaceInterpValue(
                        CubeValues, GaussCurvature, nw_pts, nw_tris, Values, i,
                        j, k, n_nw_pts, n_nw_tris);

                    Jwn += pmmc_CubeSurfaceInterpValue(
                        CubeValues, MeanCurvature, nw_pts, nw_tris, Values, i,
                        j, k, n_nw_pts, n_nw_tris);

                    Xwn += geomavg_EulerCharacteristic(
                        nw_pts, nw_tris, n_nw_pts, n_nw_tris, i, j, k);

                    // Integrate the trimmed mean curvature (hard-coded to use a distance of 4 pixels)
                    pmmc_CubeTrimSurfaceInterpValues(
                        CubeValues, MeanCurvature, SDs, nw_pts, nw_tris, Values,
                        DistanceValues, i, j, k, n_nw_pts, n_nw_tris, trimdist,
                        dummy, trJwn);

                    pmmc_CubeTrimSurfaceInterpInverseValues(
                        CubeValues, MeanCurvature, SDs, nw_pts, nw_tris, Values,
                        DistanceValues, i, j, k, n_nw_pts, n_nw_tris, trimdist,
                        dummy, trRwn);
                }
                // wns common curve averages
                if (n_local_nws_pts > 0) {
                    efawns += pmmc_CubeContactAngle(
                        CubeValues, Values, SDn_x, SDn_y, SDn_z, SDs_x, SDs_y,
                        SDs_z, local_nws_pts, i, j, k, n_local_nws_pts);

                    for (int p = 0; p < n_local_nws_pts; p++) {
                        // Extract the line segment
                        Point A = local_nws_pts(p);
                        double value = Values(p);
                        fprintf(ANGLES, "%.8g %.8g %.8g %.8g\n", A.x, A.y, A.z,
                                value);
                    }

                    pmmc_CurveCurvature(SDn, SDs, SDn_x, SDn_y, SDn_z, SDs_x,
                                        SDs_y, SDs_z, KNwns_values,
                                        KGwns_values, KNwns, KGwns, nws_pts,
                                        n_nws_pts, i, j, k);

                    lwns +=
                        pmmc_CubeCurveLength(local_nws_pts, n_local_nws_pts);

                    /* half contribution for vertices / edges at the common line 
                     *   each cube with contact line has a net of undercounting vertices 
                     *   each cube is undercounting edges due to internal counts
                    */
                    Xwn += 0.25 * n_local_nws_pts - 0.5;
                    Xws += 0.25 * n_local_nws_pts - 0.5;
                    Xns += 0.25 * n_local_nws_pts - 0.5;
                }

                // Solid interface averagees
                if (n_local_sol_tris > 0) {
                    As += pmmc_CubeSurfaceArea(local_sol_pts, local_sol_tris,
                                               n_local_sol_tris);

                    // Compute the surface orientation and the interfacial area
                    ans += pmmc_CubeSurfaceOrientation(Gns, ns_pts, ns_tris,
                                                       n_ns_tris);
                    aws += pmmc_CubeSurfaceOrientation(Gws, ws_pts, ws_tris,
                                                       n_ws_tris);

                    Xws += geomavg_EulerCharacteristic(
                        ws_pts, ws_tris, n_ws_pts, n_ws_tris, i, j, k);

                    Xns += geomavg_EulerCharacteristic(
                        ns_pts, ns_tris, n_ns_pts, n_ns_tris, i, j, k);
                }
                //...........................................................................
                // Compute the integral curvature of the non-wetting phase

                n_nw_pts = n_nw_tris = 0;
                // Compute the non-wetting phase surface and associated area
                An +=
                    geomavg_MarchingCubes(SDn, fluid_isovalue, i, j, k, nw_pts,
                                          n_nw_pts, nw_tris, n_nw_tris);
                // Compute the integral of mean curvature
                if (n_nw_pts > 0) {
                    pmmc_CubeTrimSurfaceInterpValues(
                        CubeValues, MeanCurvature, SDs, nw_pts, nw_tris, Values,
                        DistanceValues, i, j, k, n_nw_pts, n_nw_tris, trimdist,
                        trawn, dummy);
                }

                Jn += pmmc_CubeSurfaceInterpValue(CubeValues, MeanCurvature,
                                                  nw_pts, nw_tris, Values, i, j,
                                                  k, n_nw_pts, n_nw_tris);
                // Compute Euler characteristic from integral of gaussian curvature
                Kn += pmmc_CubeSurfaceInterpValue(CubeValues, GaussCurvature,
                                                  nw_pts, nw_tris, Values, i, j,
                                                  k, n_nw_pts, n_nw_tris);

                euler += geomavg_EulerCharacteristic(nw_pts, nw_tris, n_nw_pts,
                                                     n_nw_tris, i, j, k);
            }
        }
    }
    fclose(ANGLES);

    Array<char> phase_label(Nx, Ny, Nz);
    Array<double> phase_distance(Nx, Ny, Nz);
    // Analyze the wetting fluid
    for (k = 0; k < Nz; k++) {
        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                n = k * Nx * Ny + j * Nx + i;
                if (!(Dm->id[n] > 0)) {
                    // Solid phase
                    phase_label(i, j, k) = 1;
                } else if (SDn(i, j, k) < 0.0) {
                    // wetting phase
                    phase_label(i, j, k) = 0;
                } else {
                    // non-wetting phase
                    phase_label(i, j, k) = 1;
                }
                phase_distance(i, j, k) =
                    2.0 * double(phase_label(i, j, k)) - 1.0;
            }
        }
    }
    CalcDist(phase_distance, phase_label, *Dm);
    wet_morph->ComputeScalar(phase_distance, 0.f);
    //printf("generating distance at rank=%i \n",Dm->rank());
    // Analyze the wetting fluid
    for (k = 0; k < Nz; k++) {
        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                n = k * Nx * Ny + j * Nx + i;
                if (!(Dm->id[n] > 0)) {
                    // Solid phase
                    phase_label(i, j, k) = 1;
                } else if (SDn(i, j, k) < 0.0) {
                    // wetting phase
                    phase_label(i, j, k) = 1;
                } else {
                    // non-wetting phase
                    phase_label(i, j, k) = 0;
                }
                phase_distance(i, j, k) =
                    2.0 * double(phase_label(i, j, k)) - 1.0;
            }
        }
    }
    CalcDist(phase_distance, phase_label, *Dm);
    nonwet_morph->ComputeScalar(phase_distance, 0.f);
}

void TwoPhase::AssignComponentLabels() {
    //int LabelNWP=1;
    //int LabelWP=2;
    // NOTE: labeling the wetting phase components is tricky! One sandstone media had over 800,000 components
    // NumberComponents_WP = ComputeGlobalPhaseComponent(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2,Dm->rank_info,PhaseID,LabelWP,Label_WP);
    // treat all wetting phase is connected
    NumberComponents_WP = 1;
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                Label_WP(i, j, k) = 0;
                //if (SDs(i,j,k) > 0.0) PhaseID(i,j,k) = 0;
                //else if (Phase(i,j,k) > 0.0) PhaseID(i,j,k) = LabelNWP;
                //else PhaseID(i,j,k) = LabelWP;
            }
        }
    }

    // Fewer non-wetting phase features are present
    //NumberComponents_NWP = ComputeGlobalPhaseComponent(Dm->Nx-2,Dm->Ny-2,Dm->Nz-2,Dm->rank_info,PhaseID,LabelNWP,Label_NWP);
    NumberComponents_NWP = ComputeGlobalBlobIDs(
        Dm->Nx - 2, Dm->Ny - 2, Dm->Nz - 2, Dm->rank_info, SDs, SDn,
        solid_isovalue, fluid_isovalue, Label_NWP, Dm->Comm);
}

void TwoPhase::ComponentAverages() {
    int i, j, k, n;
    int kmin, kmax;
    int LabelWP, LabelNWP;
    double TempLocal;

    int cube[8][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0},
                      {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};

    ComponentAverages_WP.resize(BLOB_AVG_COUNT, NumberComponents_WP);
    ComponentAverages_NWP.resize(BLOB_AVG_COUNT, NumberComponents_NWP);

    ComponentAverages_WP.fill(0.0);
    ComponentAverages_NWP.fill(0.0);

    if (Dm->rank() == 0) {
        printf("Number of wetting phase components is %i \n",
               NumberComponents_WP);
        printf("Number of non-wetting phase components is %i \n",
               NumberComponents_NWP);
    }

    // If external boundary conditions are set, do not average over the inlet
    kmin = 1;
    kmax = Nz - 1;
    if (Dm->BoundaryCondition > 0 && Dm->kproc() == 0)
        kmin = 4;
    if (Dm->BoundaryCondition > 0 && Dm->kproc() == Dm->nprocz() - 1)
        kmax = Nz - 4;

    for (k = kmin; k < kmax; k++) {
        for (j = 1; j < Ny - 1; j++) {
            for (i = 1; i < Nx - 1; i++) {
                LabelWP = GetCubeLabel(i, j, k, Label_WP);
                LabelNWP = GetCubeLabel(i, j, k, Label_NWP);

                n_nw_pts = n_ns_pts = n_ws_pts = n_nws_pts = n_local_sol_pts =
                    n_local_nws_pts = 0;
                n_nw_tris = n_ns_tris = n_ws_tris = n_nws_seg =
                    n_local_sol_tris = 0;

                // Initialize the averaged quantities
                awn = aws = ans = lwns = 0.0;
                vawn(0) = vawn(1) = vawn(2) = 0.0;
                vawns(0) = vawns(1) = vawns(2) = 0.0;
                Gwn(0) = Gwn(1) = Gwn(2) = 0.0;
                Gwn(3) = Gwn(4) = Gwn(5) = 0.0;
                Gws(0) = Gws(1) = Gws(2) = 0.0;
                Gws(3) = Gws(4) = Gws(5) = 0.0;
                Gns(0) = Gns(1) = Gns(2) = 0.0;
                Gns(3) = Gns(4) = Gns(5) = 0.0;
                KGwns = KNwns = 0.0;
                Jwn = Kwn = efawns = 0.0;
                trawn = trJwn = 0.0;
                euler = 0.0;

                //...........................................................................
                //...........................................................................
                // Compute volume averages
                for (int p = 0; p < 8; p++) {
                    n = i + cube[p][0] + (j + cube[p][1]) * Nx +
                        (k + cube[p][2]) * Nx * Ny;
                    if (Dm->id[n] > 0) {
                        // 1-D index for this cube corner
                        // compute the norm of the gradient of the phase indicator field
                        // Compute the non-wetting phase volume contribution
                        if (Phase(i + cube[p][0], j + cube[p][1],
                                  k + cube[p][2]) > 0.0 &&
                            !(LabelNWP < 0)) {
                            // volume
                            ComponentAverages_NWP(VOL, LabelNWP) += 0.125;
                            // velocity
                            ComponentAverages_NWP(VX, LabelNWP) +=
                                0.125 * Vel_x(n);
                            ComponentAverages_NWP(VY, LabelNWP) +=
                                0.125 * Vel_y(n);
                            ComponentAverages_NWP(VZ, LabelNWP) +=
                                0.125 * Vel_z(n);
                            // center of mass
                            ComponentAverages_NWP(CMX, LabelNWP) +=
                                0.125 * (i + cube[p][0] + Dm->iproc() * Nx);
                            ComponentAverages_NWP(CMY, LabelNWP) +=
                                0.125 * (j + cube[p][1] + Dm->jproc() * Ny);
                            ComponentAverages_NWP(CMZ, LabelNWP) +=
                                0.125 * (k + cube[p][2] + Dm->kproc() * Nz);

                            // twice the kinetic energy
                            ComponentAverages_NWP(VSQ, LabelNWP) +=
                                0.125 *
                                (Vel_x(n) * Vel_x(n) + Vel_y(n) * Vel_y(n) +
                                 Vel_z(n) * Vel_z(n));

                            // volume the for pressure averaging excludes the interfacial region
                            if (DelPhi(n) < 1e-4) {
                                ComponentAverages_NWP(TRIMVOL, LabelNWP) +=
                                    0.125;
                                ComponentAverages_NWP(PRS, LabelNWP) +=
                                    0.125 * Press(n);
                            }
                        } else if (!(LabelWP < 0)) {
                            ComponentAverages_WP(VOL, LabelWP) += 0.125;
                            // velocity
                            ComponentAverages_WP(VX, LabelWP) +=
                                0.125 * Vel_x(n);
                            ComponentAverages_WP(VY, LabelWP) +=
                                0.125 * Vel_y(n);
                            ComponentAverages_WP(VZ, LabelWP) +=
                                0.125 * Vel_z(n);
                            // Center of mass
                            ComponentAverages_WP(CMX, LabelWP) +=
                                0.125 * (i + cube[p][0] + Dm->iproc() * Nx);
                            ComponentAverages_WP(CMY, LabelWP) +=
                                0.125 * (j + cube[p][1] + Dm->jproc() * Ny);
                            ComponentAverages_WP(CMZ, LabelWP) +=
                                0.125 * (k + cube[p][2] + Dm->kproc() * Nz);
                            // twice the kinetic energy
                            ComponentAverages_WP(VSQ, LabelWP) +=
                                0.125 *
                                (Vel_x(n) * Vel_x(n) + Vel_y(n) * Vel_y(n) +
                                 Vel_z(n) * Vel_z(n));

                            // volume the for pressure averaging excludes the interfacial region
                            if (DelPhi(n) < 1e-4) {
                                ComponentAverages_WP(TRIMVOL, LabelWP) += 0.125;
                                ComponentAverages_WP(PRS, LabelWP) +=
                                    0.125 * Press(n);
                            }
                        }
                    }
                }
                //...........................................................................
                // Construct the interfaces and common curve
                pmmc_ConstructLocalCube(
                    SDs, SDn, solid_isovalue, fluid_isovalue, nw_pts, nw_tris,
                    Values, ns_pts, ns_tris, ws_pts, ws_tris, local_nws_pts,
                    nws_pts, nws_seg, local_sol_pts, local_sol_tris,
                    n_local_sol_tris, n_local_sol_pts, n_nw_pts, n_nw_tris,
                    n_ws_pts, n_ws_tris, n_ns_tris, n_ns_pts, n_local_nws_pts,
                    n_nws_pts, n_nws_seg, i, j, k, Nx, Ny, Nz);

                //...........................................................................
                // wn interface averages
                if (n_nw_pts > 0 && LabelNWP >= 0 && LabelWP >= 0) {
                    // Mean curvature
                    TempLocal = pmmc_CubeSurfaceInterpValue(
                        CubeValues, MeanCurvature, nw_pts, nw_tris, Values, i,
                        j, k, n_nw_pts, n_nw_tris);
                    ComponentAverages_WP(JWN, LabelWP) += TempLocal;
                    ComponentAverages_NWP(JWN, LabelNWP) += TempLocal;

                    // Trimmed Mean curvature
                    pmmc_CubeTrimSurfaceInterpValues(
                        CubeValues, MeanCurvature, SDs, nw_pts, nw_tris, Values,
                        DistanceValues, i, j, k, n_nw_pts, n_nw_tris, trimdist,
                        trawn, trJwn);
                    ComponentAverages_WP(TRAWN, LabelWP) += trawn;
                    ComponentAverages_WP(TRJWN, LabelWP) += trJwn;
                    ComponentAverages_NWP(TRAWN, LabelNWP) += trawn;
                    ComponentAverages_NWP(TRJWN, LabelNWP) += trJwn;

                    // Gaussian curvature
                    TempLocal = pmmc_CubeSurfaceInterpValue(
                        CubeValues, GaussCurvature, nw_pts, nw_tris, Values, i,
                        j, k, n_nw_pts, n_nw_tris);
                    ComponentAverages_WP(KWN, LabelWP) += TempLocal;
                    ComponentAverages_NWP(KWN, LabelNWP) += TempLocal;

                    // Compute the normal speed of the interface
                    pmmc_InterfaceSpeed(dPdt, SDn_x, SDn_y, SDn_z, CubeValues,
                                        nw_pts, nw_tris, NormalVector,
                                        InterfaceSpeed, vawn, i, j, k, n_nw_pts,
                                        n_nw_tris);
                    ComponentAverages_WP(VWNX, LabelWP) += vawn(0);
                    ComponentAverages_WP(VWNY, LabelWP) += vawn(1);
                    ComponentAverages_WP(VWNZ, LabelWP) += vawn(2);
                    ComponentAverages_NWP(VWNX, LabelNWP) += vawn(0);
                    ComponentAverages_NWP(VWNY, LabelNWP) += vawn(1);
                    ComponentAverages_NWP(VWNZ, LabelNWP) += vawn(2);

                    // Interfacial Area
                    TempLocal = pmmc_CubeSurfaceOrientation(Gwn, nw_pts,
                                                            nw_tris, n_nw_tris);
                    ComponentAverages_WP(AWN, LabelWP) += TempLocal;
                    ComponentAverages_NWP(AWN, LabelNWP) += TempLocal;

                    ComponentAverages_WP(GWNXX, LabelWP) += Gwn(0);
                    ComponentAverages_WP(GWNYY, LabelWP) += Gwn(1);
                    ComponentAverages_WP(GWNZZ, LabelWP) += Gwn(2);
                    ComponentAverages_WP(GWNXY, LabelWP) += Gwn(3);
                    ComponentAverages_WP(GWNXZ, LabelWP) += Gwn(4);
                    ComponentAverages_WP(GWNYZ, LabelWP) += Gwn(5);

                    ComponentAverages_NWP(GWNXX, LabelNWP) += Gwn(0);
                    ComponentAverages_NWP(GWNYY, LabelNWP) += Gwn(1);
                    ComponentAverages_NWP(GWNZZ, LabelNWP) += Gwn(2);
                    ComponentAverages_NWP(GWNXY, LabelNWP) += Gwn(3);
                    ComponentAverages_NWP(GWNXZ, LabelNWP) += Gwn(4);
                    ComponentAverages_NWP(GWNYZ, LabelNWP) += Gwn(5);
                }
                //...........................................................................
                // Common curve averages
                if (n_local_nws_pts > 0 && LabelNWP >= 0 && LabelWP >= 0) {
                    // Contact angle
                    TempLocal = pmmc_CubeContactAngle(
                        CubeValues, Values, SDn_x, SDn_y, SDn_z, SDs_x, SDs_y,
                        SDs_z, local_nws_pts, i, j, k, n_local_nws_pts);
                    ComponentAverages_WP(CWNS, LabelWP) += TempLocal;
                    ComponentAverages_NWP(CWNS, LabelNWP) += TempLocal;

                    // Kinematic velocity of the common curve
                    pmmc_CommonCurveSpeed(
                        CubeValues, dPdt, vawns, SDn_x, SDn_y, SDn_z, SDs_x,
                        SDs_y, SDs_z, local_nws_pts, i, j, k, n_local_nws_pts);
                    ComponentAverages_WP(VWNSX, LabelWP) += vawns(0);
                    ComponentAverages_WP(VWNSY, LabelWP) += vawns(1);
                    ComponentAverages_WP(VWNSZ, LabelWP) += vawns(2);
                    ComponentAverages_NWP(VWNSX, LabelNWP) += vawns(0);
                    ComponentAverages_NWP(VWNSY, LabelNWP) += vawns(1);
                    ComponentAverages_NWP(VWNSZ, LabelNWP) += vawns(2);

                    // Curvature of the common curve
                    pmmc_CurveCurvature(SDn, SDs, SDn_x, SDn_y, SDn_z, SDs_x,
                                        SDs_y, SDs_z, KNwns_values,
                                        KGwns_values, KNwns, KGwns, nws_pts,
                                        n_nws_pts, i, j, k);
                    ComponentAverages_WP(KNWNS, LabelWP) += KNwns;
                    ComponentAverages_WP(KGWNS, LabelWP) += KGwns;
                    ComponentAverages_NWP(KNWNS, LabelNWP) += KNwns;
                    ComponentAverages_NWP(KGWNS, LabelNWP) += KGwns;

                    // Length of the common curve
                    TempLocal =
                        pmmc_CubeCurveLength(local_nws_pts, n_local_nws_pts);
                    ComponentAverages_NWP(LWNS, LabelNWP) += TempLocal;
                }
                //...........................................................................
                // Solid interface averages
                if (n_local_sol_pts > 0 && LabelWP >= 0) {
                    As += pmmc_CubeSurfaceArea(local_sol_pts, local_sol_tris,
                                               n_local_sol_tris);
                    // Compute the surface orientation and the interfacial area

                    TempLocal = pmmc_CubeSurfaceOrientation(Gws, ws_pts,
                                                            ws_tris, n_ws_tris);
                    ComponentAverages_WP(AWS, LabelWP) += TempLocal;
                }
                if (n_ns_pts > 0 && LabelNWP >= 0) {
                    TempLocal = pmmc_CubeSurfaceOrientation(Gns, ns_pts,
                                                            ns_tris, n_ns_tris);
                    ComponentAverages_NWP(ANS, LabelNWP) += TempLocal;
                }
                //...........................................................................

                /* Compute the Euler characteristic
				 * count all vertices, edges and faces (triangles)
				 * Euler Number = vertices - edges + faces
				 * double geomavg_EulerCharacteristic(PointList, PointCount, TriList, TriCount);
				 */
                if (LabelNWP >= 0) {
                    // find non-wetting phase interface
                    //					double euler;
                    n_nw_pts = n_nw_tris = 0;
                    geomavg_MarchingCubes(SDn, fluid_isovalue, i, j, k, nw_pts,
                                          n_nw_pts, nw_tris, n_nw_tris);
                    // Compute Euler characteristic from integral of gaussian curvature
                    euler = pmmc_CubeSurfaceInterpValue(
                        CubeValues, GaussCurvature, nw_pts, nw_tris, Values, i,
                        j, k, n_nw_pts, n_nw_tris);
                    ComponentAverages_NWP(INTCURV, LabelNWP) += euler;

                    // Compute the Euler characteristic from vertices - faces + edges
                    euler = geomavg_EulerCharacteristic(
                        nw_pts, nw_tris, n_nw_pts, n_nw_tris, i, j, k);
                    ComponentAverages_NWP(EULER, LabelNWP) += euler;
                }
            }
        }
    }

    Dm->Comm.barrier();
    if (Dm->rank() == 0) {
        printf("Component averages computed locally -- reducing result... \n");
    }
    // Globally reduce the non-wetting phase averages
    RecvBuffer.resize(BLOB_AVG_COUNT, NumberComponents_NWP);

    /*	for (int b=0; b<NumberComponents_NWP; b++){
		Dm->Comm.barrier();
		Dm->Comm.sumReduce(&ComponentAverages_NWP(0,b),&RecvBuffer(0),BLOB_AVG_COUNT);
		for (int idx=0; idx<BLOB_AVG_COUNT; idx++) ComponentAverages_NWP(idx,b)=RecvBuffer(idx);
	}
	*/
    Dm->Comm.barrier();
    Dm->Comm.sumReduce(ComponentAverages_NWP.data(), RecvBuffer.data(),
                       BLOB_AVG_COUNT * NumberComponents_NWP);
    //	Dm->Comm.sumReduce(ComponentAverages_NWP.data(),RecvBuffer.data(),BLOB_AVG_COUNT);

    if (Dm->rank() == 0) {
        printf("rescaling... \n");
    }

    for (int b = 0; b < NumberComponents_NWP; b++) {
        for (int idx = 0; idx < BLOB_AVG_COUNT; idx++)
            ComponentAverages_NWP(idx, b) = RecvBuffer(idx, b);
    }
    for (int b = 0; b < NumberComponents_NWP; b++) {
        if (ComponentAverages_NWP(VOL, b) > 0.0) {
            double Vn, pn, awn, ans, Jwn, Kwn, lwns, cwns, vsq;

            Vn = ComponentAverages_NWP(VOL, b);
            awn = ComponentAverages_NWP(AWN, b);
            ans = ComponentAverages_NWP(ANS, b);
            van(0) = ComponentAverages_NWP(VX, b) / Vn;
            van(1) = ComponentAverages_NWP(VY, b) / Vn;
            van(2) = ComponentAverages_NWP(VZ, b) / Vn;
            vsq = ComponentAverages_NWP(VSQ, b) / Vn;
            NULL_USE(ans);

            if (ComponentAverages_NWP(TRIMVOL, b) > 0.0) {
                pn = ComponentAverages_NWP(PRS, b) /
                     ComponentAverages_NWP(TRIMVOL, b);
            } else
                pn = 0.0;

            if (awn != 0.0) {
                Jwn = ComponentAverages_NWP(JWN, b) / awn;
                Kwn = ComponentAverages_NWP(KWN, b) / awn;
                vawn(0) = ComponentAverages_NWP(VWNSX, b) / awn;
                vawn(1) = ComponentAverages_NWP(VWNSY, b) / awn;
                vawn(2) = ComponentAverages_NWP(VWNSZ, b) / awn;
                Gwn(0) = ComponentAverages_NWP(GWNXX, b) / awn;
                Gwn(1) = ComponentAverages_NWP(GWNYY, b) / awn;
                Gwn(2) = ComponentAverages_NWP(GWNZZ, b) / awn;
                Gwn(3) = ComponentAverages_NWP(GWNXY, b) / awn;
                Gwn(4) = ComponentAverages_NWP(GWNXZ, b) / awn;
                Gwn(5) = ComponentAverages_NWP(GWNYZ, b) / awn;
            } else
                Jwn = Kwn = 0.0;

            trawn = ComponentAverages_NWP(TRAWN, b);
            trJwn = ComponentAverages_NWP(TRJWN, b);
            if (trawn > 0.0)
                trJwn /= trawn;

            lwns = ComponentAverages_NWP(LWNS, b);
            if (lwns != 0.0) {
                cwns = ComponentAverages_NWP(CWNS, b) / lwns;
                vawns(0) = ComponentAverages_NWP(VWNSX, b) / lwns;
                vawns(1) = ComponentAverages_NWP(VWNSY, b) / lwns;
                vawns(2) = ComponentAverages_NWP(VWNSZ, b) / lwns;
            } else
                cwns = 0.0;

            ComponentAverages_NWP(PRS, b) = pn;
            ComponentAverages_NWP(VX, b) = van(0);
            ComponentAverages_NWP(VY, b) = van(1);
            ComponentAverages_NWP(VZ, b) = van(2);
            ComponentAverages_NWP(VSQ, b) = vsq;

            ComponentAverages_NWP(JWN, b) = Jwn;
            ComponentAverages_NWP(KWN, b) = Kwn;
            ComponentAverages_NWP(VWNX, b) = vawn(0);
            ComponentAverages_NWP(VWNY, b) = vawn(1);
            ComponentAverages_NWP(VWNZ, b) = vawn(2);

            ComponentAverages_NWP(GWNXX, b) = Gwn(0);
            ComponentAverages_NWP(GWNYY, b) = Gwn(1);
            ComponentAverages_NWP(GWNZZ, b) = Gwn(2);
            ComponentAverages_NWP(GWNXY, b) = Gwn(3);
            ComponentAverages_NWP(GWNXZ, b) = Gwn(4);
            ComponentAverages_NWP(GWNYZ, b) = Gwn(5);

            ComponentAverages_NWP(CWNS, b) = cwns;
            ComponentAverages_NWP(VWNSX, b) = vawns(0);
            ComponentAverages_NWP(VWNSY, b) = vawns(1);
            ComponentAverages_NWP(VWNSZ, b) = vawns(2);

            ComponentAverages_NWP(CMX, b) /= Vn;
            ComponentAverages_NWP(CMY, b) /= Vn;
            ComponentAverages_NWP(CMZ, b) /= Vn;

            ComponentAverages_NWP(TRJWN, b) = trJwn;

            ComponentAverages_NWP(EULER, b) /= (2 * PI);
        }
    }

    if (Dm->rank() == 0) {
        printf("reduce WP averages... \n");
    }

    // reduce the wetting phase averages
    for (int b = 0; b < NumberComponents_WP; b++) {
        Dm->Comm.barrier();
        Dm->Comm.sumReduce(&ComponentAverages_WP(0, b), RecvBuffer.data(),
                           BLOB_AVG_COUNT);
        for (int idx = 0; idx < BLOB_AVG_COUNT; idx++)
            ComponentAverages_WP(idx, b) = RecvBuffer(idx);
    }

    for (int b = 0; b < NumberComponents_WP; b++) {
        if (ComponentAverages_WP(VOL, b) > 0.0) {
            double Vw, pw, awn, ans, Jwn, Kwn, lwns, cwns, vsq;
            Vw = ComponentAverages_WP(VOL, b);
            awn = ComponentAverages_WP(AWN, b);
            ans = ComponentAverages_WP(ANS, b);
            vaw(0) = ComponentAverages_WP(VX, b) / Vw;
            vaw(1) = ComponentAverages_WP(VY, b) / Vw;
            vaw(2) = ComponentAverages_WP(VZ, b) / Vw;
            vsq = ComponentAverages_WP(VSQ, b) / Vw;
            NULL_USE(ans);

            if (ComponentAverages_WP(TRIMVOL, b) > 0.0) {
                pw = ComponentAverages_WP(PRS, b) /
                     ComponentAverages_WP(TRIMVOL, b);
            } else
                pw = 0.0;

            if (awn != 0.0) {
                Jwn = ComponentAverages_WP(JWN, b) / awn;
                Kwn = ComponentAverages_WP(KWN, b) / awn;
                vawn(0) = ComponentAverages_WP(VWNSX, b) / awn;
                vawn(1) = ComponentAverages_WP(VWNSY, b) / awn;
                vawn(2) = ComponentAverages_WP(VWNSZ, b) / awn;
                Gwn(0) = ComponentAverages_WP(GWNXX, b) / awn;
                Gwn(1) = ComponentAverages_WP(GWNYY, b) / awn;
                Gwn(2) = ComponentAverages_WP(GWNZZ, b) / awn;
                Gwn(3) = ComponentAverages_WP(GWNXY, b) / awn;
                Gwn(4) = ComponentAverages_WP(GWNXZ, b) / awn;
                Gwn(5) = ComponentAverages_WP(GWNYZ, b) / awn;
            } else
                Jwn = Kwn = 0.0;

            trawn = ComponentAverages_WP(TRAWN, b);
            trJwn = ComponentAverages_WP(TRJWN, b);
            if (trawn > 0.0)
                trJwn /= trawn;

            lwns = ComponentAverages_WP(LWNS, b);
            if (lwns != 0.0) {
                cwns = ComponentAverages_WP(CWNS, b) / lwns;
                vawns(0) = ComponentAverages_WP(VWNSX, b) / lwns;
                vawns(1) = ComponentAverages_WP(VWNSY, b) / lwns;
                vawns(2) = ComponentAverages_WP(VWNSZ, b) / lwns;
            } else
                cwns = 0.0;

            ComponentAverages_WP(PRS, b) = pw;
            ComponentAverages_WP(VX, b) = vaw(0);
            ComponentAverages_WP(VY, b) = vaw(1);
            ComponentAverages_WP(VZ, b) = vaw(2);
            ComponentAverages_WP(VSQ, b) = vsq;

            ComponentAverages_WP(JWN, b) = Jwn;
            ComponentAverages_WP(KWN, b) = Kwn;
            ComponentAverages_WP(VWNX, b) = vawn(0);
            ComponentAverages_WP(VWNY, b) = vawn(1);
            ComponentAverages_WP(VWNZ, b) = vawn(2);

            ComponentAverages_WP(GWNXX, b) = Gwn(0);
            ComponentAverages_WP(GWNYY, b) = Gwn(1);
            ComponentAverages_WP(GWNZZ, b) = Gwn(2);
            ComponentAverages_WP(GWNXY, b) = Gwn(3);
            ComponentAverages_WP(GWNXZ, b) = Gwn(4);
            ComponentAverages_WP(GWNYZ, b) = Gwn(5);

            ComponentAverages_WP(CWNS, b) = cwns;
            ComponentAverages_WP(VWNSX, b) = vawns(0);
            ComponentAverages_WP(VWNSY, b) = vawns(1);
            ComponentAverages_WP(VWNSZ, b) = vawns(2);

            ComponentAverages_WP(TRJWN, b) = trJwn;
        }
    }
}

void TwoPhase::Reduce() {
    int i;
    //double iVol_global = 1.0 / Volume;
    //...........................................................................
    Dm->Comm.barrier();
    nwp_volume_global = Dm->Comm.sumReduce(nwp_volume);
    wp_volume_global = Dm->Comm.sumReduce(wp_volume);
    awn_global = Dm->Comm.sumReduce(awn);
    ans_global = Dm->Comm.sumReduce(ans);
    aws_global = Dm->Comm.sumReduce(aws);
    lwns_global = Dm->Comm.sumReduce(lwns);
    As_global = Dm->Comm.sumReduce(As);
    Jwn_global = Dm->Comm.sumReduce(Jwn);
    Kwn_global = Dm->Comm.sumReduce(Kwn);
    KGwns_global = Dm->Comm.sumReduce(KGwns);
    KNwns_global = Dm->Comm.sumReduce(KNwns);
    efawns_global = Dm->Comm.sumReduce(efawns);
    wwndnw_global = Dm->Comm.sumReduce(wwndnw);
    wwnsdnwn_global = Dm->Comm.sumReduce(wwnsdnwn);
    Jwnwwndnw_global = Dm->Comm.sumReduce(Jwnwwndnw);
    // Phase averages
    vol_w_global = Dm->Comm.sumReduce(vol_w);
    vol_n_global = Dm->Comm.sumReduce(vol_n);
    paw_global = Dm->Comm.sumReduce(paw);
    pan_global = Dm->Comm.sumReduce(pan);
    for (int idx = 0; idx < 3; idx++)
        vaw_global(idx) = Dm->Comm.sumReduce(vaw(idx));
    for (int idx = 0; idx < 3; idx++)
        van_global(idx) = Dm->Comm.sumReduce(van(idx));
    for (int idx = 0; idx < 3; idx++)
        vawn_global(idx) = Dm->Comm.sumReduce(vawn(idx));
    for (int idx = 0; idx < 3; idx++)
        vawns_global(idx) = Dm->Comm.sumReduce(vawns(idx));
    for (int idx = 0; idx < 6; idx++) {
        Gwn_global(idx) = Dm->Comm.sumReduce(Gwn(idx));
        Gns_global(idx) = Dm->Comm.sumReduce(Gns(idx));
        Gws_global(idx) = Dm->Comm.sumReduce(Gws(idx));
    }
    trawn_global = Dm->Comm.sumReduce(trawn);
    trJwn_global = Dm->Comm.sumReduce(trJwn);
    trRwn_global = Dm->Comm.sumReduce(trRwn);
    Xwn_global = Dm->Comm.sumReduce(Xwn);
    Xws_global = Dm->Comm.sumReduce(Xws);
    Xns_global = Dm->Comm.sumReduce(Xns);
    An_global = Dm->Comm.sumReduce(An);
    Jn_global = Dm->Comm.sumReduce(Jn);
    Kn_global = Dm->Comm.sumReduce(Kn);
    euler_global = Dm->Comm.sumReduce(euler);

    Dm->Comm.barrier();

    // Normalize the phase averages
    // (density of both components = 1.0)
    if (vol_w_global > 0.0) {
        paw_global = paw_global / vol_w_global;
    }
    if (wp_volume_global > 0.0) {
        vaw_global(0) = vaw_global(0) / wp_volume_global;
        vaw_global(1) = vaw_global(1) / wp_volume_global;
        vaw_global(2) = vaw_global(2) / wp_volume_global;
    }
    if (vol_n_global > 0.0) {
        pan_global = pan_global / vol_n_global;
    }
    if (nwp_volume_global > 0.0) {
        van_global(0) = van_global(0) / nwp_volume_global;
        van_global(1) = van_global(1) / nwp_volume_global;
        van_global(2) = van_global(2) / nwp_volume_global;
    }
    // Normalize surface averages by the interfacial area
    if (awn_global > 0.0) {
        Jwn_global /= awn_global;
        //Kwn_global /= awn_global;
        wwndnw_global /= awn_global;
        for (i = 0; i < 3; i++)
            vawn_global(i) /= awn_global;
        for (i = 0; i < 6; i++)
            Gwn_global(i) /= awn_global;
    }

    if (lwns_global > 0.0) {
        efawns_global /= lwns_global;
        //KNwns_global /= lwns_global;
        //KGwns_global /= lwns_global;
        for (i = 0; i < 3; i++)
            vawns_global(i) /= lwns_global;
    }
    if (trawn_global > 0.0) {
        trJwn_global /= trawn_global;
        trRwn_global /= trawn_global;
        trRwn_global = 2.0 * fabs(trRwn_global);
        trJwn_global = fabs(trJwn_global);
    }

    if (ans_global > 0.0)
        for (i = 0; i < 6; i++)
            Gns_global(i) /= ans_global;
    if (aws_global > 0.0)
        for (i = 0; i < 6; i++)
            Gws_global(i) /= aws_global;

    euler_global /= (2 * PI);
    sat_w = 1.0 - nwp_volume_global / (nwp_volume_global + wp_volume_global);

    // Compute the specific interfacial areas and common line length (dimensionless per unit volume)
    /* 
    awn_global = awn_global * iVol_global;
    ans_global = ans_global * iVol_global;
    aws_global = aws_global * iVol_global;
    dEs = dEs * iVol_global;
    lwns_global = lwns_global * iVol_global;
    */
}

void TwoPhase::NonDimensionalize(double D, double viscosity, double IFT) {
    NULL_USE(viscosity);
    NULL_USE(IFT);
    awn_global *= D;
    ans_global *= D;
    ans_global *= D;
    lwns_global *= D * D;
}

void TwoPhase::PrintStatic() {
    if (Dm->rank() == 0) {
        FILE *STATIC;
        STATIC = fopen("geometry.csv", "a+");
        if (fseek(STATIC, 0, SEEK_SET) == fseek(STATIC, 0, SEEK_CUR)) {
            // If timelog is empty, write a short header to list the averages
            fprintf(STATIC, "sw awn ans aws Jwn Kwn lwns cwns KGws "
                            "KGwn Xwn Xws Xns "); // Scalar averages
            fprintf(
                STATIC,
                "Gwnxx Gwnyy Gwnzz Gwnxy Gwnxz Gwnyz "); // Orientation tensors
            fprintf(STATIC, "Gwsxx Gwsyy Gwszz Gwsxy Gwsxz Gwsyz ");
            fprintf(STATIC, "Gnsxx Gnsyy Gnszz Gnsxy Gnsxz Gnsyz ");
            fprintf(STATIC, "trawn trJwn trRwn "); //trimmed curvature,
            fprintf(STATIC, "Vw Aw Jw Xw ");       //miknowski measures,
            fprintf(STATIC, "Vn An Jn Xn\n");      //miknowski measures,
            //fprintf(STATIC,"Euler Kn2 Jn2 An2\n"); 			//miknowski measures,
        }

        fprintf(STATIC, "%.5g ", sat_w); // saturation
        fprintf(STATIC, "%.5g %.5g %.5g ", awn_global, ans_global,
                aws_global); // interfacial areas
        fprintf(STATIC, "%.5g %.5g ", Jwn_global,
                Kwn_global);                     // curvature of wn interface
        fprintf(STATIC, "%.5g ", lwns_global);   // common curve length
        fprintf(STATIC, "%.5g ", efawns_global); // average contact angle
        fprintf(STATIC, "%.5g %.5g ", KNwns_global,
                KGwns_global); // total curvature contributions of common line
        fprintf(STATIC, "%.5g %.5g %.5g ", Xwn_global, Xns_global,
                Xws_global); // Euler characteristic
        fprintf(STATIC, "%.5g %.5g %.5g %.5g %.5g %.5g ", Gwn_global(0),
                Gwn_global(1), Gwn_global(2), Gwn_global(3), Gwn_global(4),
                Gwn_global(5)); // orientation of wn interface
        fprintf(STATIC, "%.5g %.5g %.5g %.5g %.5g %.5g ", Gns_global(0),
                Gns_global(1), Gns_global(2), Gns_global(3), Gns_global(4),
                Gns_global(5)); // orientation of ns interface
        fprintf(STATIC, "%.5g %.5g %.5g %.5g %.5g %.5g ", Gws_global(0),
                Gws_global(1), Gws_global(2), Gws_global(3), Gws_global(4),
                Gws_global(5)); // orientation of ws interface
        fprintf(STATIC, "%.5g %.5g %.5g ", trawn_global, trJwn_global,
                trRwn_global); // Trimmed curvature
        fprintf(STATIC, "%.5g %.5g %.5g %.5g ", wet_morph->V(), wet_morph->A(),
                wet_morph->H(), wet_morph->X());
        fprintf(STATIC, "%.5g %.5g %.5g %.5g\n", nonwet_morph->V(),
                nonwet_morph->A(), nonwet_morph->H(), nonwet_morph->X());
        //fprintf(STATIC,"%.5g %.5g %.5g %.5g\n",euler_global, Kn_global, Jn_global, An_global);			// minkowski measures
        fclose(STATIC);
    }
}

void TwoPhase::PrintAll(int timestep) {
    if (Dm->rank() == 0) {
        fprintf(TIMELOG, "%i %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g ",
                timestep, rho_n, rho_w, nu_n, nu_w, Fx, Fy, Fz, gamma_wn);
        fprintf(TIMELOG, "%.5g %.5g %.5g ", sat_w, paw_global,
                pan_global); // saturation and pressure
        fprintf(TIMELOG, "%.5g %.5g %.5g ", awn_global, ans_global,
                aws_global); // interfacial areas
        fprintf(TIMELOG, "%.5g %.5g ", Jwn_global,
                Kwn_global);                      // curvature of wn interface
        fprintf(TIMELOG, "%.5g ", lwns_global);   // common curve length
        fprintf(TIMELOG, "%.5g ", efawns_global); // average contact angle
        fprintf(TIMELOG, "%.5g %.5g ", KNwns_global,
                KGwns_global); // curvature of wn interface
        fprintf(TIMELOG, "%.5g %.5g %.5g ", vaw_global(0), vaw_global(1),
                vaw_global(2)); // average velocity of w phase
        fprintf(TIMELOG, "%.5g %.5g %.5g ", van_global(0), van_global(1),
                van_global(2)); // average velocity of n phase
        fprintf(TIMELOG, "%.5g %.5g %.5g ", vawn_global(0), vawn_global(1),
                vawn_global(2)); // velocity of wn interface
        fprintf(TIMELOG, "%.5g %.5g %.5g ", vawns_global(0), vawns_global(1),
                vawns_global(2)); // velocity of wn interface
        fprintf(TIMELOG, "%.5g %.5g %.5g %.5g %.5g %.5g ", Gwn_global(0),
                Gwn_global(1), Gwn_global(2), Gwn_global(3), Gwn_global(4),
                Gwn_global(5)); // orientation of wn interface
        fprintf(TIMELOG, "%.5g %.5g %.5g %.5g %.5g %.5g ", Gns_global(0),
                Gns_global(1), Gns_global(2), Gns_global(3), Gns_global(4),
                Gns_global(5)); // orientation of ns interface
        fprintf(TIMELOG, "%.5g %.5g %.5g %.5g %.5g %.5g ", Gws_global(0),
                Gws_global(1), Gws_global(2), Gws_global(3), Gws_global(4),
                Gws_global(5)); // orientation of ws interface
        fprintf(TIMELOG, "%.5g %.5g %.5g ", trawn_global, trJwn_global,
                trRwn_global); // Trimmed curvature
        fprintf(TIMELOG, "%.5g %.5g %.5g ", wwndnw_global, wwnsdnwn_global,
                Jwnwwndnw_global); // kinematic quantities
        fprintf(TIMELOG, "%.5g %.5g %.5g %.5g ", wet_morph->V(), wet_morph->A(),
                wet_morph->H(), wet_morph->X());
        fprintf(TIMELOG, "%.5g %.5g %.5g %.5g\n", nonwet_morph->V(),
                nonwet_morph->A(), nonwet_morph->H(), nonwet_morph->X());
        //		fprintf(TIMELOG,"%.5g %.5g %.5g %.5g\n",euler_global, Kn_global, Jn_global, An_global);			// minkowski measures
        fflush(TIMELOG);
    } else {
        sat_w = 1.0 - nwp_volume / (nwp_volume + wp_volume);

        fprintf(TIMELOG, "%i %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g ",
                timestep, rho_n, rho_w, nu_n, nu_w, Fx, Fy, Fz, gamma_wn);
        fprintf(TIMELOG, "%.5g %.5g %.5g ", sat_w, paw,
                pan); // saturation and pressure
        fprintf(TIMELOG, "%.5g %.5g %.5g ", awn, ans, aws); // interfacial areas
        fprintf(TIMELOG, "%.5g %.5g ", Jwn, Kwn); // curvature of wn interface
        fprintf(TIMELOG, "%.5g ", lwns);          // common curve length
        fprintf(TIMELOG, "%.5g ", efawns);        // average contact angle
        fprintf(TIMELOG, "%.5g %.5g ", KNwns,
                KGwns); // curvature of wn interface
        fprintf(TIMELOG, "%.5g %.5g %.5g ", vaw(0), vaw(1),
                vaw(2)); // average velocity of w phase
        fprintf(TIMELOG, "%.5g %.5g %.5g ", van(0), van(1),
                van(2)); // average velocity of n phase
        fprintf(TIMELOG, "%.5g %.5g %.5g ", vawn(0), vawn(1),
                vawn(2)); // velocity of wn interface
        fprintf(TIMELOG, "%.5g %.5g %.5g ", vawns(0), vawns(1),
                vawns(2)); // velocity of wn interface
        fprintf(TIMELOG, "%.5g %.5g %.5g %.5g %.5g %.5g ", Gwn(0), Gwn(1),
                Gwn(2), Gwn(3), Gwn(4), Gwn(5)); // orientation of wn interface
        fprintf(TIMELOG, "%.5g %.5g %.5g %.5g %.5g %.5g ", Gns(0), Gns(1),
                Gns(2), Gns(3), Gns(4), Gns(5)); // orientation of ns interface
        fprintf(TIMELOG, "%.5g %.5g %.5g %.5g %.5g %.5g ", Gws(0), Gws(1),
                Gws(2), Gws(3), Gws(4), Gws(5)); // orientation of ws interface
        fprintf(TIMELOG, "%.5g %.5g %.5g ", trawn, trJwn,
                trRwn); // Trimmed curvature
        fprintf(TIMELOG, "%.5g %.5g %.5g ", wwndnw, wwnsdnwn,
                Jwnwwndnw); // kinematic quantities
        fprintf(TIMELOG, "%.5g %.5g %.5g %.5g ", wet_morph->Vi, wet_morph->Ai,
                wet_morph->Ji, wet_morph->Xi);
        fprintf(TIMELOG, "%.5g %.5g %.5g %.5g\n", nonwet_morph->Vi,
                nonwet_morph->Ai, nonwet_morph->Ji, nonwet_morph->Xi);
        //		fprintf(TIMELOG,"%.5g %.5g %.5g %.5g\n",euler, Kn, Jn, An);			// minkowski measures
        fflush(TIMELOG);
    }
}

void TwoPhase::PrintComponents(int timestep) {
    if (Dm->rank() == 0) {
        printf("PRINT %i COMPONENT AVEREAGES: time = %i \n",
               (int)ComponentAverages_NWP.size(1), timestep);
        for (int b = 0; b < NumberComponents_NWP; b++) {
            //if (ComponentAverages_NWP(TRIMVOL,b) > 0.0){
            fprintf(NWPLOG, "%i ", timestep - 5);
            if (Label_NWP_map.empty()) {
                // print index
                fprintf(NWPLOG, "%i ", b);
            } else {
                // print final global id
                fprintf(NWPLOG, "%i ", Label_NWP_map[b]);
            }
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(VOL, b));
            //			fprintf(NWPLOG,"%.5g ",ComponentAverages_NWP(TRIMVOL,b));
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(PRS, b));
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(AWN, b));
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(ANS, b));
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(JWN, b));
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(KWN, b));
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(LWNS, b));
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(CWNS, b));
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(VX, b));
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(VY, b));
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(VZ, b));
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(VWNX, b));
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(VWNY, b));
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(VWNZ, b));
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(VWNSX, b));
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(VWNSY, b));
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(VWNSZ, b));
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(VSQ, b));
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(GWNXX, b));
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(GWNYY, b));
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(GWNZZ, b));
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(GWNXY, b));
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(GWNXZ, b));
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(GWNYZ, b));
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(CMX, b));
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(CMY, b));
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(CMZ, b));
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(TRAWN, b));
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(TRJWN, b));
            fprintf(NWPLOG, "%.5g ", ComponentAverages_NWP(INTCURV, b));
            fprintf(NWPLOG, "%.5g\n", ComponentAverages_NWP(EULER, b));
            //				fprintf(NWPLOG,"%.5g ",ComponentAverages_NWP(NVERT,b));
            //			fprintf(NWPLOG,"%.5g ",ComponentAverages_NWP(NSIDE,b));
            //		fprintf(NWPLOG,"%.5g\n",ComponentAverages_NWP(NFACE,b));
            //}
        }
        fflush(NWPLOG);

        for (int b = 0; b < NumberComponents_WP; b++) {
            if (ComponentAverages_WP(TRIMVOL, b) > 0.0) {
                fprintf(WPLOG, "%i ", timestep - 5);
                fprintf(WPLOG, "%i ", b);
                fprintf(WPLOG, "%.5g ", ComponentAverages_WP(VOL, b));
                //			fprintf(WPLOG,"%.5g ",ComponentAverages_WP(TRIMVOL,b));
                fprintf(WPLOG, "%.5g ", ComponentAverages_WP(PRS, b));
                fprintf(WPLOG, "%.5g ", ComponentAverages_WP(AWN, b));
                fprintf(WPLOG, "%.5g ", ComponentAverages_WP(AWS, b));
                fprintf(WPLOG, "%.5g ", ComponentAverages_WP(JWN, b));
                fprintf(WPLOG, "%.5g ", ComponentAverages_WP(KWN, b));
                fprintf(WPLOG, "%.5g ", ComponentAverages_WP(LWNS, b));
                fprintf(WPLOG, "%.5g ", ComponentAverages_WP(CWNS, b));
                fprintf(WPLOG, "%.5g ", ComponentAverages_WP(VX, b));
                fprintf(WPLOG, "%.5g ", ComponentAverages_WP(VY, b));
                fprintf(WPLOG, "%.5g ", ComponentAverages_WP(VZ, b));
                fprintf(WPLOG, "%.5g ", ComponentAverages_WP(VWNX, b));
                fprintf(WPLOG, "%.5g ", ComponentAverages_WP(VWNY, b));
                fprintf(WPLOG, "%.5g ", ComponentAverages_WP(VWNZ, b));
                fprintf(WPLOG, "%.5g ", ComponentAverages_WP(VWNSX, b));
                fprintf(WPLOG, "%.5g ", ComponentAverages_WP(VWNSY, b));
                fprintf(WPLOG, "%.5g ", ComponentAverages_WP(VWNSZ, b));
                fprintf(WPLOG, "%.5g ", ComponentAverages_WP(VSQ, b));
                fprintf(WPLOG, "%.5g ", ComponentAverages_WP(GWNXX, b));
                fprintf(WPLOG, "%.5g ", ComponentAverages_WP(GWNYY, b));
                fprintf(WPLOG, "%.5g ", ComponentAverages_WP(GWNZZ, b));
                fprintf(WPLOG, "%.5g ", ComponentAverages_WP(GWNXY, b));
                fprintf(WPLOG, "%.5g ", ComponentAverages_WP(GWNXZ, b));
                fprintf(WPLOG, "%.5g ", ComponentAverages_WP(GWNYZ, b));
                fprintf(WPLOG, "%.5g ", ComponentAverages_WP(CMX, b));
                fprintf(WPLOG, "%.5g ", ComponentAverages_WP(CMY, b));
                fprintf(WPLOG, "%.5g ", ComponentAverages_WP(CMZ, b));
                fprintf(WPLOG, "%.5g ", ComponentAverages_WP(TRAWN, b));
                fprintf(WPLOG, "%.5g\n", ComponentAverages_WP(TRJWN, b));
            }
        }
        fflush(WPLOG);
    }
}

inline int TwoPhase::GetCubeLabel(int i, int j, int k, IntArray &BlobLabel) {
    int label;
    label = BlobLabel(i, j, k);
    label = max(label, BlobLabel(i + 1, j, k));
    label = max(label, BlobLabel(i, j + 1, k));
    label = max(label, BlobLabel(i + 1, j + 1, k));
    label = max(label, BlobLabel(i, j, k + 1));
    label = max(label, BlobLabel(i + 1, j, k + 1));
    label = max(label, BlobLabel(i, j + 1, k + 1));
    label = max(label, BlobLabel(i + 1, j + 1, k + 1));

    return label;
}

void TwoPhase::SortBlobs() {
    //printf("Sorting the blobs based on volume \n");
    //printf("-----------------------------------------------\n");
    int TempLabel, a, aa, bb, i, j, k, idx;
    double TempValue;
    //.......................................................................
    // Sort NWP components by volume
    //.......................................................................
    IntArray OldLabel(NumberComponents_NWP);
    for (a = 0; a < NumberComponents_NWP; a++)
        OldLabel(a) = a;
    // Sort the blob averages based on volume
    for (aa = 0; aa < NumberComponents_NWP - 1; aa++) {
        for (bb = aa + 1; bb < NumberComponents_NWP; bb++) {
            if (ComponentAverages_NWP(VOL, aa) <
                ComponentAverages_NWP(VOL, bb)) {
                // Exchange location of blobs aa and bb
                //printf("Switch blob %i with %i \n", OldLabel(aa),OldLabel(bb));
                // switch the label
                TempLabel = OldLabel(bb);
                OldLabel(bb) = OldLabel(aa);
                OldLabel(aa) = TempLabel;
                // switch the averages
                for (idx = 0; idx < BLOB_AVG_COUNT; idx++) {
                    TempValue = ComponentAverages_NWP(idx, bb);
                    ComponentAverages_NWP(idx, bb) =
                        ComponentAverages_NWP(idx, aa);
                    ComponentAverages_NWP(idx, aa) = TempValue;
                }
            }
        }
    }
    IntArray NewLabel(NumberComponents_NWP);
    for (aa = 0; aa < NumberComponents_NWP; aa++) {
        // Match the new label for original blob aa
        bb = 0;
        while (OldLabel(bb) != aa)
            bb++;
        NewLabel(aa) = bb;
    }
    // Re-label the blob ID
    //	printf("Re-labeling the blobs, now indexed by volume \n");
    for (k = 0; k < Nz; k++) {
        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                if (Label_NWP(i, j, k) > -1) {
                    TempLabel = NewLabel(Label_NWP(i, j, k));
                    Label_NWP(i, j, k) = TempLabel;
                }
            }
        }
    }
    //.......................................................................
    // Sort WP components by volume
    //.......................................................................
    OldLabel.resize(NumberComponents_WP);
    for (a = 0; a < NumberComponents_WP; a++)
        OldLabel(a) = a;
    // Sort the blob averages based on volume
    for (aa = 0; aa < NumberComponents_WP - 1; aa++) {
        for (bb = aa + 1; bb < NumberComponents_WP; bb++) {
            if (ComponentAverages_WP(VOL, aa) < ComponentAverages_WP(VOL, bb)) {
                // Exchange location of blobs aa and bb
                //printf("Switch blob %i with %i \n", OldLabel(aa),OldLabel(bb));
                // switch the label
                TempLabel = OldLabel(bb);
                OldLabel(bb) = OldLabel(aa);
                OldLabel(aa) = TempLabel;
                // switch the averages
                for (idx = 0; idx < BLOB_AVG_COUNT; idx++) {
                    TempValue = ComponentAverages_WP(idx, bb);
                    ComponentAverages_WP(idx, bb) =
                        ComponentAverages_WP(idx, aa);
                    ComponentAverages_WP(idx, aa) = TempValue;
                }
            }
        }
    }
    NewLabel.resize(NumberComponents_WP);
    for (aa = 0; aa < NumberComponents_WP; aa++) {
        // Match the new label for original blob aa
        bb = 0;
        while (OldLabel(bb) != aa)
            bb++;
        NewLabel(aa) = bb;
    }
    // Re-label the blob ID
    //	printf("Re-labeling the blobs, now indexed by volume \n");
    for (k = 0; k < Nz; k++) {
        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                if (Label_WP(i, j, k) > -1) {
                    TempLabel = NewLabel(Label_WP(i, j, k));
                    Label_WP(i, j, k) = TempLabel;
                }
            }
        }
    }
    //.......................................................................
}

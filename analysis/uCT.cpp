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
#include "analysis/uCT.h"
#include "analysis/analysis.h"
#include "analysis/distance.h"
#include "analysis/filters.h"
#include "analysis/imfilter.h"

template <class T> inline int sign(T x) {
    if (x == 0)
        return 0;
    return x > 0 ? 1 : -1;
}

inline float trilinear(float dx, float dy, float dz, float f1, float f2,
                       float f3, float f4, float f5, float f6, float f7,
                       float f8) {
    double f, dx2, dy2, dz2, h0, h1;
    dx2 = 1.0 - dx;
    dy2 = 1.0 - dy;
    dz2 = 1.0 - dz;
    h0 = (dx * f2 + dx2 * f1) * dy2 + (dx * f4 + dx2 * f3) * dy;
    h1 = (dx * f6 + dx2 * f5) * dy2 + (dx * f8 + dx2 * f7) * dy;
    f = h0 * dz2 + h1 * dz;
    return (f);
}

void InterpolateMesh(const Array<float> &Coarse, Array<float> &Fine) {
    PROFILE_START("InterpolateMesh");

    // Interpolate values from a Coarse mesh to a fine one
    // This routine assumes cell-centered meshes with 1 ghost cell

    // Fine mesh
    int Nx = int(Fine.size(0)) - 2;
    int Ny = int(Fine.size(1)) - 2;
    int Nz = int(Fine.size(2)) - 2;

    // Coarse mesh
    int nx = int(Coarse.size(0)) - 2;
    int ny = int(Coarse.size(1)) - 2;
    int nz = int(Coarse.size(2)) - 2;

    // compute the stride
    int hx = Nx / nx;
    int hy = Ny / ny;
    int hz = Nz / nz;
    ASSERT(nx * hx == Nx);
    ASSERT(ny * hy == Ny);
    ASSERT(nz * hz == Nz);

    // value to map distance between meshes (since distance is in voxels)
    //  usually hx=hy=hz (or something very close)
    //  the mapping is not exact
    //  however, it's assumed the coarse solution will be refined
    //  a good guess is the goal here!
    float mapvalue = sqrt(hx * hx + hy * hy + hz * hz);

    // Interpolate to the fine mesh
    for (int k = -1; k < Nz + 1; k++) {
        int k0 = floor((k - 0.5 * hz) / hz);
        int k1 = k0 + 1;
        int k2 = k0 + 2;
        float dz = ((k + 0.5) - (k0 + 0.5) * hz) / hz;
        ASSERT(k0 >= -1 && k0 < nz + 1 && dz >= 0 && dz <= 1);
        for (int j = -1; j < Ny + 1; j++) {
            int j0 = floor((j - 0.5 * hy) / hy);
            int j1 = j0 + 1;
            int j2 = j0 + 2;
            float dy = ((j + 0.5) - (j0 + 0.5) * hy) / hy;
            ASSERT(j0 >= -1 && j0 < ny + 1 && dy >= 0 && dy <= 1);
            for (int i = -1; i < Nx + 1; i++) {
                int i0 = floor((i - 0.5 * hx) / hx);
                int i1 = i0 + 1;
                int i2 = i0 + 2;
                float dx = ((i + 0.5) - (i0 + 0.5) * hx) / hx;
                ASSERT(i0 >= -1 && i0 < nx + 1 && dx >= 0 && dx <= 1);
                float val = trilinear(
                    dx, dy, dz, Coarse(i1, j1, k1), Coarse(i2, j1, k1),
                    Coarse(i1, j2, k1), Coarse(i2, j2, k1), Coarse(i1, j1, k2),
                    Coarse(i2, j1, k2), Coarse(i1, j2, k2), Coarse(i2, j2, k2));
                Fine(i + 1, j + 1, k + 1) = mapvalue * val;
            }
        }
    }
    PROFILE_STOP("InterpolateMesh");
}

// Smooth the data using the distance
void smooth(const Array<float> &VOL, const Array<float> &Dist, float sigma,
            Array<float> &MultiScaleSmooth, fillHalo<float> &fillFloat) {
    for (size_t i = 0; i < VOL.length(); i++) {
        // use exponential weight based on the distance
        float dst = Dist(i);
        float tmp = exp(-(dst * dst) / (sigma * sigma));
        float value = dst > 0 ? -1 : 1;
        MultiScaleSmooth(i) = tmp * VOL(i) + (1 - tmp) * value;
    }
    fillFloat.fill(MultiScaleSmooth);
}

// Segment the data
void segment(const Array<float> &data, Array<char> &ID, float tol) {
    ASSERT(data.size() == ID.size());
    for (size_t i = 0; i < data.length(); i++) {
        if (data(i) > tol)
            ID(i) = 0;
        else
            ID(i) = 1;
    }
}

// Remove disconnected phases
void removeDisconnected(Array<char> &ID, const Domain &Dm) {
    // Run blob identification to remove disconnected volumes
    BlobIDArray GlobalBlobID;
    DoubleArray SignDist(ID.size());
    DoubleArray Phase(ID.size());
    for (size_t i = 0; i < ID.length(); i++) {
        SignDist(i) = (2 * ID(i) - 1);
        Phase(i) = 1;
    }
    ComputeGlobalBlobIDs(ID.size(0) - 2, ID.size(1) - 2, ID.size(2) - 2,
                         Dm.rank_info, Phase, SignDist, 0, 0, GlobalBlobID,
                         Dm.Comm);
    for (size_t i = 0; i < ID.length(); i++) {
        if (GlobalBlobID(i) > 0)
            ID(i) = 0;
        ID(i) = GlobalBlobID(i);
    }
}

// Solve a level (without any coarse level information)
void solve(const Array<float> &VOL, Array<float> &Mean, Array<char> &ID,
           Array<float> &Dist, Array<float> &MultiScaleSmooth,
           Array<float> &NonLocalMean, fillHalo<float> &fillFloat,
           const Domain &Dm, int nprocx, float threshold, float lamda,
           float sigsq, int depth) {
    PROFILE_SCOPED(timer, "solve");
    // Compute the median filter on the sparse array
    Med3D(VOL, Mean);
    fillFloat.fill(Mean);
    segment(Mean, ID, threshold);
    // Compute the distance using the segmented volume
    CalcDist(Dist, ID, Dm);
    fillFloat.fill(Dist);
    smooth(VOL, Dist, 2.0, MultiScaleSmooth, fillFloat);
    // Compute non-local mean
    //	int depth = 5;
    //	float sigsq=0.1;
    int nlm_count =
        NLM3D(MultiScaleSmooth, Mean, Dist, NonLocalMean, depth, sigsq);
    NULL_USE(nlm_count);
    fillFloat.fill(NonLocalMean);
}

// Refine a solution from a coarse grid to a fine grid
void refine(const Array<float> &Dist_coarse, const Array<float> &VOL,
            Array<float> &Mean, Array<char> &ID, Array<float> &Dist,
            Array<float> &MultiScaleSmooth, Array<float> &NonLocalMean,
            fillHalo<float> &fillFloat, const Domain &Dm, int nprocx, int level,
            float threshold, float lamda, float sigsq, int depth) {
    PROFILE_SCOPED(timer, "refine");
    int ratio[3] = {int(Dist.size(0) / Dist_coarse.size(0)),
                    int(Dist.size(1) / Dist_coarse.size(1)),
                    int(Dist.size(2) / Dist_coarse.size(2))};
    // Interpolate the distance from the coarse to fine grid
    InterpolateMesh(Dist_coarse, Dist);
    // Compute the median filter on the array and segment
    Med3D(VOL, Mean);
    fillFloat.fill(Mean);
    segment(Mean, ID, threshold);
    // If the ID has the wrong distance, set the distance to 0 and run a simple filter to set neighbors to 0
    for (size_t i = 0; i < ID.length(); i++) {
        char id = Dist(i) > 0 ? 1 : 0;
        if (id != ID(i))
            Dist(i) = 0;
    }
    fillFloat.fill(Dist);
    std::function<float(int, const float *)> filter_1D = [](int N,
                                                            const float *data) {
        bool zero = data[0] == 0 || data[2] == 0;
        return zero ? data[1] * 1e-12 : data[1];
    };
    std::vector<imfilter::BC> BC(3, imfilter::BC::replicate);
    std::vector<std::function<float(int, const float *)>> filter_set(3,
                                                                     filter_1D);
    Dist = imfilter::imfilter_separable<float>(Dist, {1, 1, 1}, filter_set, BC);
    fillFloat.fill(Dist);
    // Smooth the volume data
    float h = 2 * lamda *
              sqrt(double(ratio[0] * ratio[0] + ratio[1] * ratio[1] +
                          ratio[2] * ratio[2]));
    smooth(VOL, Dist, h, MultiScaleSmooth, fillFloat);
    // Compute non-local mean
    //	int depth = 3;
    //	float sigsq = 0.1;
    int nlm_count =
        NLM3D(MultiScaleSmooth, Mean, Dist, NonLocalMean, depth, sigsq);
    NULL_USE(nlm_count);
    fillFloat.fill(NonLocalMean);
    segment(NonLocalMean, ID, 0.001);
    for (size_t i = 0; i < ID.length(); i++) {
        char id = Dist(i) > 0 ? 1 : 0;
        if (id != ID(i) || fabs(Dist(i)) < 1)
            Dist(i) = 2.0 * ID(i) - 1.0;
    }
    // Remove disconnected domains
    //removeDisconnected( ID, Dm );
    // Compute the distance using the segmented volume
    if (level > 0) {
        CalcDist(Dist, ID, Dm);
        fillFloat.fill(Dist);
    }
}

// Remove regions that are likely noise by shrinking the volumes by dx,
// removing all values that are more than dx+delta from the surface, and then
// growing by dx+delta and intersecting with the original data
void filter_final(Array<char> &ID, Array<float> &Dist,
                  fillHalo<float> &fillFloat, const Domain &Dm,
                  Array<float> &Mean, Array<float> &Dist1,
                  Array<float> &Dist2) {
    PROFILE_SCOPED(timer, "filter_final");
    int rank = Dm.Comm.getRank();
    int Nx = Dm.Nx - 2;
    int Ny = Dm.Ny - 2;
    int Nz = Dm.Nz - 2;
    // Calculate the distance
    CalcDist(Dist, ID, Dm);
    fillFloat.fill(Dist);
    // Compute the range to shrink the volume based on the L2 norm of the distance
    Array<float> Dist0(Nx, Ny, Nz);
    fillFloat.copy(Dist, Dist0);
    float tmp = 0;
    for (size_t i = 0; i < Dist0.length(); i++)
        tmp += Dist0(i) * Dist0(i);
    tmp =
        sqrt(Dm.Comm.sumReduce(tmp) / Dm.Comm.sumReduce<float>(Dist0.length()));
    const float dx1 = 0.3 * tmp;
    const float dx2 = 1.05 * dx1;
    if (rank == 0)
        printf("   %0.1f %0.1f %0.1f\n", tmp, dx1, dx2);
    // Update the IDs/Distance removing regions that are < dx of the range
    Dist1 = Dist;
    Dist2 = Dist;
    Array<char> ID1 = ID;
    Array<char> ID2 = ID;
    for (size_t i = 0; i < ID.length(); i++) {
        ID1(i) = Dist(i) < -dx1 ? 1 : 0;
        ID2(i) = Dist(i) > dx1 ? 1 : 0;
    }
    //Array<float> Dist1 = Dist;
    //Array<float> Dist2 = Dist;
    CalcDist(Dist1, ID1, Dm);
    CalcDist(Dist2, ID2, Dm);
    fillFloat.fill(Dist1);
    fillFloat.fill(Dist2);
    // Keep those regions that are within dx2 of the new volumes
    Mean = Dist;
    for (size_t i = 0; i < ID.length(); i++) {
        if (Dist1(i) + dx2 > 0 && ID(i) <= 0) {
            Mean(i) = -1;
        } else if (Dist2(i) + dx2 > 0 && ID(i) > 0) {
            Mean(i) = 1;
        } else {
            Mean(i) = Dist(i) > 0 ? 0.5 : -0.5;
        }
    }
    // Find regions of uncertainty that are entirely contained within another region
    fillHalo<double> fillDouble(Dm.Comm, Dm.rank_info, {Nx, Ny, Nz}, {1, 1, 1},
                                0, 1);
    fillHalo<BlobIDType> fillInt(Dm.Comm, Dm.rank_info, {Nx, Ny, Nz}, {1, 1, 1},
                                 0, 1);
    BlobIDArray GlobalBlobID;
    DoubleArray SignDist(ID.size());
    for (size_t i = 0; i < ID.length(); i++)
        SignDist(i) = fabs(Mean(i)) == 1 ? -1 : 1;
    fillDouble.fill(SignDist);
    DoubleArray Phase(ID.size());
    Phase.fill(1);
    ComputeGlobalBlobIDs(Nx, Ny, Nz, Dm.rank_info, Phase, SignDist, 0, 0,
                         GlobalBlobID, Dm.Comm);
    fillInt.fill(GlobalBlobID);
    int N_blobs = Dm.Comm.maxReduce(GlobalBlobID.max() + 1);
    std::vector<float> mean(N_blobs, 0);
    std::vector<int> count(N_blobs, 0);
    for (int k = 1; k <= Nz; k++) {
        for (int j = 1; j <= Ny; j++) {
            for (int i = 1; i <= Nx; i++) {
                int id = GlobalBlobID(i, j, k);
                if (id >= 0) {
                    if (GlobalBlobID(i - 1, j, k) < 0) {
                        mean[id] += Mean(i - 1, j, k);
                        count[id]++;
                    }
                    if (GlobalBlobID(i + 1, j, k) < 0) {
                        mean[id] += Mean(i + 1, j, k);
                        count[id]++;
                    }
                    if (GlobalBlobID(i, j - 1, k) < 0) {
                        mean[id] += Mean(i, j - 1, k);
                        count[id]++;
                    }
                    if (GlobalBlobID(i, j + 1, k) < 0) {
                        mean[id] += Mean(i, j + 1, k);
                        count[id]++;
                    }
                    if (GlobalBlobID(i, j, k - 1) < 0) {
                        mean[id] += Mean(i, j, k - 1);
                        count[id]++;
                    }
                    if (GlobalBlobID(i, j, k + 1) < 0) {
                        mean[id] += Mean(i, j, k + 1);
                        count[id]++;
                    }
                }
            }
        }
    }
    mean = Dm.Comm.sumReduce(mean);
    count = Dm.Comm.sumReduce(count);
    for (size_t i = 0; i < mean.size(); i++)
        mean[i] /= count[i];
    /*if (rank==0) {
        for (size_t i=0; i<mean.size(); i++)
            printf("%i %0.4f\n",i,mean[i]);
    }*/
    for (size_t i = 0; i < Mean.length(); i++) {
        int id = GlobalBlobID(i);
        if (id >= 0) {
            if (fabs(mean[id]) > 0.95) {
                // Isolated domain surrounded by one domain
                GlobalBlobID(i) = -2;
                Mean(i) = sign(mean[id]);
            } else {
                // Boarder volume, set to liquid
                Mean(i) = 1;
            }
        }
    }
    // Perform the final segmentation and update the distance
    fillFloat.fill(Mean);
    segment(Mean, ID, 0.01);
    CalcDist(Dist, ID, Dm);
    fillFloat.fill(Dist);
}

// Filter the original data
void filter_src(const Domain &Dm, Array<float> &src) {
    PROFILE_START("Filter source data");
    int Nx = Dm.Nx - 2;
    int Ny = Dm.Ny - 2;
    int Nz = Dm.Nz - 2;
    fillHalo<float> fillFloat(Dm.Comm, Dm.rank_info, {Nx, Ny, Nz}, {1, 1, 1}, 0,
                              1);
    // Perform a hot-spot filter on the data
    std::vector<imfilter::BC> BC = {imfilter::BC::replicate,
                                    imfilter::BC::replicate,
                                    imfilter::BC::replicate};
    std::function<float(const Array<float> &)> filter_3D =
        [](const Array<float> &data) {
            float min1 = std::min(data(0, 1, 1), data(2, 1, 1));
            float min2 = std::min(data(1, 0, 1), data(1, 2, 1));
            float min3 = std::min(data(1, 1, 0), data(1, 1, 2));
            float max1 = std::max(data(0, 1, 1), data(2, 1, 1));
            float max2 = std::max(data(1, 0, 1), data(1, 2, 1));
            float max3 = std::max(data(1, 1, 0), data(1, 1, 2));
            float min = std::min(min1, std::min(min2, min3));
            float max = std::max(max1, std::max(max2, max3));
            return std::max(std::min(data(1, 1, 1), max), min);
        };
    std::function<float(const Array<float> &)> filter_1D =
        [](const Array<float> &data) {
            float min = std::min(data(0), data(2));
            float max = std::max(data(0), data(2));
            return std::max(std::min(data(1), max), min);
        };
    //LOCVOL[0] = imfilter::imfilter<float>( LOCVOL[0], {1,1,1}, filter_3D, BC );
    std::vector<std::function<float(const Array<float> &)>> filter_set(
        3, filter_1D);
    src = imfilter::imfilter_separable<float>(src, {1, 1, 1}, filter_set, BC);
    fillFloat.fill(src);
    // Perform a gaussian filter on the data
    int Nh[3] = {2, 2, 2};
    float sigma[3] = {1.0, 1.0, 1.0};
    std::vector<Array<float>> H(3);
    H[0] = imfilter::create_filter<float>({Nh[0]}, "gaussian", &sigma[0]);
    H[1] = imfilter::create_filter<float>({Nh[1]}, "gaussian", &sigma[1]);
    H[2] = imfilter::create_filter<float>({Nh[2]}, "gaussian", &sigma[2]);
    src = imfilter::imfilter_separable(src, H, BC);
    fillFloat.fill(src);
    PROFILE_STOP("Filter source data");
}

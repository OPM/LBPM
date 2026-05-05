/*
<<<<<<< HEAD
  Copyright 2026
  Diogo Nardelli Siebert, Universidade Federal de Santa Catarina
  Bernardo Gehlen, Universidade Federal de Santa Catarina
  Alexandre Miers Zabot, Universidade Federal de Santa Catarina

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the :contentReference[oaicite:0]{index=0}, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM. If not, see <http://www.gnu.org/licenses/>.
*/

/*
  ----------------------------------------------------------------------
  Methodological Notes

  This file implements the method described in:

  A. M. Zabot et al.,
  "A Unified Algorithm for the Young–Laplace Method Applied to Porous Media,"
  Brazilian Journal of Physics, vol. 54, no. 3, 2024, p. 63.

  The distance transform is computed using the method from:

  P. F. Felzenszwalb and D. P. Huttenlocher,
  "Distance transforms of sampled functions,"
  Theory of Computing, vol. 8, no. 1, pp. 415–428, 2012.

  Connected component labeling is based on:

  L. He, Y. Chao, and K. Suzuki,
  "A run-based two-scan labeling algorithm,"
  IEEE Transactions on Image Processing, vol. 17, no. 5, pp. 749–756, 2008.
  ----------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <stdint.h>
#include <string>
#include <math.h>
#include <bits/stdc++.h>
#include "../common/Array.h"
#include "../common/Domain.h"
#include "../common/UtilityMacros.h"

#define SOLID ((unsigned char)0)
#define DISPLACED ((unsigned char)1)
#define INJECTED ((unsigned char)2)

#define BACKGROUND 0
#define FOREGROUND 1

using namespace std;

void merge(const int &u, const int &v, vector<int> &next, vector<int> &tail,
           vector<int> &rtable) {
    for (int i = v; i != -1;) {
        rtable[i] = u;
        i = next[i];
    }
    next[tail[u]] = v;
    tail[u] = tail[v];
}

void resolve(const int &x, const int &y, vector<int> &next, vector<int> &tail,
             vector<int> &rtable) {
    const int u = rtable[x];
    const int v = rtable[y];
    if (u < v)
        merge(u, v, next, tail, rtable);
    else if (v < u)
        merge(v, u, next, tail, rtable);
}

void component_labeling(IntArray &IMG, const int &F, const int &B) {

    size_t maxNumberOfLabels = (IMG.length() + 1) / 2;

    std::vector<int> next(maxNumberOfLabels);
    std::vector<int> tail(maxNumberOfLabels);
    std::vector<int> rtable(maxNumberOfLabels);

    const int nx = IMG.size(0);
    const int ny = IMG.size(1);
    const int nz = IMG.size(2);

    int lx = 0, nl = 1;
    vector<int> uniq_labels(3);
    int nuniq;

    for (int z = 0; z < nz; z++) {
        for (int y = 0; y < ny; y++) {
            for (int x = 0; x < nx; x++) {

                if (IMG(x, y, z) == F) {
                    const int lq = (x > 0) ? IMG(x - 1, y, z) : B;
                    const int lp = (y > 0) ? IMG(x, y - 1, z) : B;
                    const int lz = (z > 0) ? IMG(x, y, z - 1) : B;

                    nuniq = 0;
                    if (lp != B) {
                        uniq_labels[nuniq] = lp;
                        nuniq++;
                    }
                    if (lq != B && lq != lp) {
                        uniq_labels[nuniq] = lq;
                        nuniq++;
                    }
                    if (lz != B && lz != lp && lz != lq) {
                        uniq_labels[nuniq] = lz;
                        nuniq++;
                    }

                    // Deals with non unique labels in neighboaring points by unifying them
                    switch (nuniq) {

                    case 0:
                        nl++;
                        lx = nl;

                        rtable[nl] = nl;
                        next[nl] = -1;
                        tail[nl] = nl;
                        break;

                    case 1:
                        lx = uniq_labels[0];
                        break;

                    case 2:

                        resolve(uniq_labels[0], uniq_labels[1], next, tail,
                                rtable);

                        lx = min(uniq_labels[0], uniq_labels[1]);
                        break;

                    case 3:

                        resolve(uniq_labels[0], uniq_labels[1], next, tail,
                                rtable);
                        resolve(uniq_labels[0], uniq_labels[2], next, tail,
                                rtable);
                        resolve(uniq_labels[1], uniq_labels[2], next, tail,
                                rtable);

                        lx = min(uniq_labels[0],
                                 min(uniq_labels[1], uniq_labels[2]));
                        break;
                    }

                    IMG(x, y, z) = lx;
                }
            }
        }
    }

    int *img = IMG.data();

    for (size_t n = 0; n < IMG.length(); n++) {
        if (*img != B)
            *img = rtable[*img];
        img++;
    }
}

// Compute the intersection between two lower-envelope parabolas.
static inline float intersection(int q, int vk, float fq, float fvk) {
    float qq = (float)q * (float)q;
    float vv = (float)vk * (float)vk;
    return ((fq + qq) - (fvk + vv)) / (2.0 * ((float)q - (float)vk));
}

// Compute the exact squared Euclidean distance transform for a 1D line.
void edt_1d(const int *f, int *g, int n) {
    int *v = (int *)malloc((size_t)n * sizeof(int));
    float *z = (float *)malloc((size_t)(n + 1) * sizeof(float));

    int k = 0;
    v[0] = 0;
    z[0] = -INFINITY;
    z[1] = INFINITY;

    for (int q = 1; q < n; q++) {
        float s;

        while (1) {
            int vk = v[k];
            s = intersection(q, vk, f[q], f[vk]);

            if (s <= z[k]) {
                k--;
                if (k < 0) {
                    k = 0;
                    break;
                }
            } else {
                break;
            }
        }

        if (k == 0) {
            int vk = v[k];
            s = intersection(q, vk, f[q], f[vk]);

            if (s <= z[k]) {
                v[0] = q;
                z[0] = -INFINITY;
                z[1] = INFINITY;
                continue;
            }
        }

        k++;
        v[k] = q;
        z[k] = s;
        z[k + 1] = INFINITY;
    }

    k = 0;
    for (int x = 0; x < n; x++) {
        while (z[k + 1] < (float)x) {
            k++;
        }

        int vk = v[k];
        int dx = x - vk;
        g[x] = dx * dx + f[vk];
    }

    free(v);
    free(z);
}

// Copy one strided line, transform it in 1D, and write it back.
void process_line(int *edt2, const size_t first, const size_t stride, int n) {

    IntArray f(n);
    IntArray g(n);

    size_t pos = first;
    for (int i = 0; i < n; i++) {
        f(i) = edt2[pos];
        pos += stride;
    }

    edt_1d(f.data(), g.data(), n);

    pos = first;
    for (int i = 0; i < n; i++) {
        edt2[pos] = g(i);
        pos += stride;
    }
}

template <typename TYPE>
void edt_3d(unsigned char target, Array<TYPE> &image, IntArray &distance2) {
    const int nx = image.size(0);
    const int ny = image.size(1);
    const int nz = image.size(2);

    size_t nvox = image.length();
    int BIG = static_cast<float>(nx * nx + ny * ny + nz * nz) + 1;

    TYPE *img = image.data();
    int *edt2 = distance2.data();

    // Initialize target voxels with zero distance and all others with a large value.
    for (int i = 0; i < (int)nvox; i++) {
        edt2[i] = (img[i] == target) ? 0 : BIG;
    }

    for (int z = 0; z < nz; z++) {
        for (int y = 0; y < ny; y++) {
            process_line(edt2, nx * (y + z * ny), 1, nx);
        }
    }

    for (int z = 0; z < nz; z++) {
        for (int x = 0; x < nx; x++) {
            process_line(edt2, z * nx * ny + x, nx, ny);
        }
    }

    for (int y = 0; y < ny; y++) {
        for (int x = 0; x < nx; x++) {
            process_line(edt2, y * nx + x, nx * ny, nz);
        }
    }
}

// Validate database options and report the accepted values on failure.
void checkOption(std::string a, std::vector<string> s, std::string keyName) {
    std::string message = "Error: Invalid option '" + a + "' for " + keyName +
                          ". Valid options are: ";
    for (string b : s) {
        if (a == b)
            return;
        message += "'" + b + "', ";
    }
    message.pop_back();
    ERROR(message);
}

template <typename TYPE>
void setRegion(Array<TYPE> &A, TYPE value, int x0, int x1, int y0, int y1,
               int z0, int z1) {
    for (int z = z0; z < z1; z++)
        for (int y = y0; y < y1; y++)
            for (int x = x0; x < x1; x++) {
                A(x, y, z) = value;
            }
}

class Full_Morphology {
public:
    Full_Morphology(int, char *[]);
    int calc(const int &);

public:
    std::vector<int> r_ini = {0, 0, 0};
    std::vector<int> r_end = {0, 0, 0};

    std::vector<int> flowAxis = {
        false, false, false}; // Injection occurs along the selected axis
    bool flowPos = false;     // Sense of invasion

    int ny, nx, nz;          // Work dimensions (with reservoirs, if any)
    int dimy, dimx, dimz;    // Original dimensions
    int NP;                  // Porous pixels
    int inletPos, outletPos; // Reservoir regions

    int rChamberI[3];
    int rChamberO[3];

    double resolution;

    vector<int> diameter; // Invasion diameters

    bool compressible = false; // Compressibility, enabled for MICP

    bool allFaces; // Add reservoirs around all active faces for MICP
    bool saveImg;

    Array<u_int8_t> originalState; // Original image, used for reference
    Array<u_int8_t> currentState;  // Work image, modified during the simulation
    Array<int16_t> finalMap;
    IntArray originalEDT;
    Array<bool> trapped;
};

Full_Morphology::Full_Morphology(int argc, char *argv[]) {

    string filename;

    filename = argv[1];

    auto db = std::make_shared<Database>(filename);

    auto domain_db = db->getDatabase("Domain");
    auto fm_db = db->getDatabase("FM");

    auto size = domain_db->getVector<int>("N");
    nx = size[0];
    ny = size[1];
    nz = size[2];

    finalMap.resize(size[0], size[1], size[2]);
    finalMap.fill(-1);

    auto ReadValues = domain_db->getVector<int>("ReadValues");
    auto WriteValues = domain_db->getVector<int>("WriteValues");

    resolution = domain_db->getScalar<double>("voxel_length");

    auto READFILE = domain_db->getScalar<std::string>("Filename");
    const string mmfile(READFILE);

    saveImg = fm_db->getScalar<bool>("SaveImage");

    auto protocol = fm_db->getScalar<std::string>("protocol");
    checkOption(protocol, {"micp", "drainage"}, "protocol");

    if (protocol == "micp") {
        compressible = true;
        allFaces = true;
        flowPos = true;
        if (size[2] > 1)
            flowAxis[2] = true;
        else if (size[1] > 1)
            flowAxis[1] = true;
        else if (size[0] > 1)
            flowAxis[0] = true;
    } else if (protocol == "drainage") {
	allFaces = false;    
        auto direction = fm_db->getScalar<std::string>("direction");
        checkOption(direction, {"+x", "-x", "+y", "-y", "+z", "-z"},
                    "direction");
        flowAxis[0] = (direction[1] == 'x');
        flowAxis[1] = (direction[1] == 'y');
        flowAxis[2] = (direction[1] == 'z');
        flowPos = (direction[0] == '+');
        compressible = fm_db->getWithDefault<bool>("compressible", false);
    }

    auto diameterRange = fm_db->getVector<int>("Diameters");
    int numberOfDiameters =
        (diameterRange[1] - diameterRange[0]) / diameterRange[2] + 1;

    if (numberOfDiameters <= 0) {
        ERROR("Error: It was impossible to create diameters array. ");
    }

    diameter.resize(numberOfDiameters);

    for (int i = 0; i < numberOfDiameters; i++) {
        diameter[i] = diameterRange[1] - i * diameterRange[2];
    }

    r_end = size;
    dimx = nx;
    dimy = ny;
    dimz = nz;

    // Add additional layers for input/output reservoirs.
    for (int i = 0; i < 3; i++) {
        if (((flowAxis[i]) && (!allFaces)) || ((size[i] > 1) && (allFaces))) {
            size[i] += 2;
            r_ini[i] = 1;
            r_end[i] = size[i] - 1;
        }
    }

    nx = size[0];
    ny = size[1];
    nz = size[2];

    originalState.resize(nx, ny, nz);
    currentState.resize(originalState.size());
    originalEDT.resize(originalState.size());
    trapped.resize(originalState.size());

    currentState.fill(INJECTED);
    trapped.fill(false);

    // For directional injection, set the outlet face to DISPLACED fluid.
    for (int i = 0; i < 3; i++) {
        int rMin[3] = {0, 0, 0};
        int rMax[3] = {nx, ny, nz};

        if (!allFaces) {
            if (flowAxis[i]) {
                inletPos = flowPos ? 0 : size[i] - 1;
                outletPos = flowPos ? size[i] - 1 : 0;
                rMin[i] = outletPos;
                rMax[i] = rMin[i] + 1;
                setRegion(currentState, DISPLACED, rMin[0], rMax[0], rMin[1],
                          rMax[1], rMin[2], rMax[2]);
            }
   	    rChamberI[i] = flowAxis[i] ? inletPos : rMax[i] / 2;
            rChamberO[i] = flowAxis[i] ? outletPos : rMax[i] / 2;
        }

        else
            rChamberI[i] = 0;
    }

    int mapValue[255] = {-1};
    for (size_t idx = 0; idx < ReadValues.size(); idx++) {

        if ((ReadValues[idx] < 0) || (ReadValues[idx] > 255)) {
            ERROR("Only values between 0 - 255 can be used as labels in "
                  "ReadValues");
            cout << ReadValues[idx] << endl;
        }
        if ((WriteValues[idx] < 0) || (WriteValues[idx] > 2)) {
            ERROR("Only values between 0 (SOLID), 1 and 2 (INJECT/DISPLACED "
                  "FLUIDS) can be used as labels in WriteValues");
        }
        mapValue[ReadValues[idx]] = (int)WriteValues[idx];
    }

    FILE *rawFile = fopen(mmfile.c_str(), "r");
    if (rawFile == NULL) {
        ERROR("Error openning the file " + mmfile);
    }

    long SEEK_BEGIN = ftell(rawFile);
    long expectedSize = (long)(r_end[2] - r_ini[2]) *
                        (long)(r_end[1] - r_ini[1]) *
                        (long)(r_end[0] - r_ini[0]);

    fseek(rawFile, 0, SEEK_END);

    if (ftell(rawFile) != expectedSize) {
        ERROR("File '" + mmfile + "' size is different from the expected (" +
              to_string(expectedSize) + " bytes).");
    }

    fseek(rawFile, 0,
          SEEK_BEGIN); // Move to the beginning of the file before reading.

    unsigned char readValue;

    NP = 0;
    for (int z = r_ini[2]; z < r_end[2]; z++) {
        for (int y = r_ini[1]; y < r_end[1]; y++) {
            for (int x = r_ini[0]; x < r_end[0]; x++) {

                fread(&readValue, sizeof(unsigned char), 1, rawFile);
                if (mapValue[readValue] == -1) {
                    ERROR(std::string("Not specified value in '" + filename +
                                      "' at (" + to_string(x) + ", " +
                                      to_string(y) + ", " + to_string(z) +
                                      ")."));
                }

                currentState(x, y, z) = (unsigned char)mapValue[readValue];
                if (currentState(x, y, z) != SOLID)
                    NP++;
            }
        }
    }

    fclose(rawFile);
    originalState = currentState;

    // Calculate the distance transform of the original solid phase.
    edt_3d(SOLID, originalState, originalEDT);

    // Create the output CSV file and header when needed.
    bool WriteHeader = false;
    FILE *log_file = fopen("injection_output.csv", "r");
    if (log_file != NULL)
        fclose(log_file);
    else
        WriteHeader = true;

    if (WriteHeader) {
        log_file = fopen("injection_output.csv", "a+");
        fprintf(log_file, "step diameter_px diameter_um num_px_in frac_in "
                          "num_px_out frac_out\n");
        fclose(log_file);
    }
}

int Full_Morphology::calc(const int &step) {

    const int D = diameter[step];
    const double D24 = D * D / 4.0;

    IntArray auxMatrix(nx, ny, nz);

    // Mark the candidate invaded region as the region where the center of a
    // D diameter sphere can be placed (erosion of the solid region) united with the reservoir.
    for (int z = 0; z < nz; z++) {
        for (int y = 0; y < ny; y++) {
            for (int x = 0; x < nx; x++) { 
 		auxMatrix(x, y, z) = (originalEDT(x, y, z) >= D24) || (originalState(x, y, z) == INJECTED) ? FOREGROUND : BACKGROUND;
            }
        }
    }

    component_labeling(auxMatrix, FOREGROUND, BACKGROUND);

    // Extract the label associated with the injection layer.
    int chamber_label = auxMatrix(rChamberI[0], rChamberI[1], rChamberI[2]);

    // Filter only the voxels connected to the input where the center of a a sphere with the specified diameter
    // can be placed.
    for (int z = 0; z < nz; z++) {
        for (int y = 0; y < ny; y++) {
            for (int x = 0; x < nx; x++) {
		    auxMatrix(x, y, z) = (originalEDT(x, y, z) >= D24 && auxMatrix(x,y,z) == chamber_label) ? FOREGROUND : BACKGROUND;
            }
        }
    }

    // Compute the distance from the filtered region
    edt_3d(FOREGROUND, auxMatrix, auxMatrix);

    // Perform the dilation of the filtered region and performs the union of this result
    // with the previous step stored at the currentState
    for (int z = 0; z < nz; z++) {
        for (int y = 0; y < ny; y++) {
            for (int x = 0; x < nx; x++) {
                if (auxMatrix(x, y, z) < D24) {
                    currentState(x, y, z) = INJECTED;
                }
      		auxMatrix(x, y, z) = (currentState(x, y, z) == INJECTED)
                                         ? FOREGROUND
                                         : BACKGROUND;
            }
        }
    }


    component_labeling(auxMatrix, FOREGROUND, BACKGROUND);
    chamber_label = auxMatrix(rChamberI[0], rChamberI[1], rChamberI[2]);

    // Readds the DISPLACED fluid output layer
    if (!allFaces) {

        if (flowAxis[0]) {
            setRegion(currentState, DISPLACED, outletPos, outletPos + 1, 0, ny,
                      0, nz);
        } else if (flowAxis[1]) {
            setRegion(currentState, DISPLACED, 0, nx, outletPos, outletPos + 1,
                      0, nz);
        } else if (flowAxis[2]) {
            setRegion(currentState, DISPLACED, 0, nx, 0, ny, outletPos,
                      outletPos + 1);
        }
    }

    // For incompressible flow, disconnected displaced regions remain trapped
    // in the final state.
    if (!compressible) {

        for (int x = 0; x < nx; x++) {
            for (int y = 0; y < ny; y++) {
                for (int z = 0; z < nz; z++) {

                    if (trapped(x, y, z))
                        currentState(x, y, z) = DISPLACED;

                    auxMatrix(x, y, z) = (currentState(x, y, z) == DISPLACED)
                                             ? FOREGROUND
                                             : BACKGROUND;
                }
            }
        }

        component_labeling(auxMatrix, FOREGROUND, BACKGROUND);
        chamber_label = auxMatrix(rChamberO[0], rChamberO[1], rChamberO[2]);

        for (int x = 0; x < nx; x++) {
            for (int y = 0; y < ny; y++) {
                for (int z = 0; z < nz; z++) {
                    if (currentState(x, y, z) == DISPLACED &&
                        auxMatrix(x, y, z) != chamber_label)
                        trapped(x, y, z) = true;
                }
            }
        }
    }

    int injectedVolume = 0, displacedVolume = 0;
    for (int z = r_ini[2]; z < r_end[2]; z++) {
        for (int y = r_ini[1]; y < r_end[1]; y++) {
            for (int x = r_ini[0]; x < r_end[0]; x++) {

                if (currentState(x, y, z) == INJECTED)
                    injectedVolume++;
                else if (currentState(x, y, z) == DISPLACED)
                    displacedVolume++;

                int16_t *value =
                    &finalMap(x - r_ini[0], y - r_ini[1], z - r_ini[2]);
                if (*value == -1 && currentState(x, y, z) == INJECTED)
                    *value = (int16_t)D;
                if (*value == -1 && originalState(x, y, z) == SOLID)
                    *value = 0;
            }
        }
    }

    if (saveImg && (D == diameter.back())) {
        FILE *FRAW;

        FRAW = fopen("invasion_diameters.raw", "wb");
        fwrite(finalMap.data(), sizeof(int16_t), finalMap.length(), FRAW);
        fclose(FRAW);

        FILE *FMHD = fopen("invasion_diameters.mhd", "w");

        fprintf(FMHD, "ObjectType = Image\n");
        fprintf(FMHD, "NDims = 3\n");
        fprintf(FMHD, "DimSize = %d %d %d\n", dimx, dimy, dimz);
        fprintf(FMHD, "ElementType =  MET_SHORT\n");
        fprintf(FMHD, "ElementSpacing = %.1f %.1f %.1f\n", resolution,
                resolution, resolution);
        fprintf(FMHD, "ElementByteOrderMSB = False\n");
        fprintf(FMHD, "ElementDataFile = %s\n", "invasion_diameters.raw");
        fprintf(FMHD, "HeaderSize = 0\n");
        fclose(FMHD);
    }

    FILE *log_file = fopen("injection_output.csv", "a");
    fprintf(log_file, "%d %d %f %d %f %d %f\n", step, D, D * resolution,
            injectedVolume, injectedVolume / (1.0 * NP), displacedVolume,
            displacedVolume / (1.0 * NP));
    fclose(log_file);

    return D;
}

int main(int argc, char *argv[]) {

    if (argc != 2)
        ERROR("Wrong number of parameters.");

    Full_Morphology fm(argc, argv);

    const int nsteps = fm.diameter.size();
    for (int step = 0; step < nsteps; step++) {

        int d = fm.calc(step);
        cout << "Step " << step << ", D = " << d << " px." << endl;
    }

    return 0;
}

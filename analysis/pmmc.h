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
#ifndef pmmc_INC
#define pmmc_INC

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "common/Array.h"
#include "PointList.h"
#include "common/Utilities.h"

using namespace std;

static int edgeTable[256] = {
    0x0,   0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c, 0x80c, 0x905, 0xa0f,
    0xb06, 0xc0a, 0xd03, 0xe09, 0xf00, 0x190, 0x99,  0x393, 0x29a, 0x596, 0x49f,
    0x795, 0x69c, 0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90, 0x230,
    0x339, 0x33,  0x13a, 0x636, 0x73f, 0x435, 0x53c, 0xa3c, 0xb35, 0x83f, 0x936,
    0xe3a, 0xf33, 0xc39, 0xd30, 0x3a0, 0x2a9, 0x1a3, 0xaa,  0x7a6, 0x6af, 0x5a5,
    0x4ac, 0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0, 0x460, 0x569,
    0x663, 0x76a, 0x66,  0x16f, 0x265, 0x36c, 0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a,
    0x963, 0xa69, 0xb60, 0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff,  0x3f5, 0x2fc,
    0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0, 0x650, 0x759, 0x453,
    0x55a, 0x256, 0x35f, 0x55,  0x15c, 0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53,
    0x859, 0x950, 0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc,  0xfcc,
    0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0, 0x8c0, 0x9c9, 0xac3, 0xbca,
    0xcc6, 0xdcf, 0xec5, 0xfcc, 0xcc,  0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9,
    0x7c0, 0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c, 0x15c, 0x55,
    0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650, 0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6,
    0xfff, 0xcf5, 0xdfc, 0x2fc, 0x3f5, 0xff,  0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c, 0x36c, 0x265, 0x16f,
    0x66,  0x76a, 0x663, 0x569, 0x460, 0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af,
    0xaa5, 0xbac, 0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa,  0x1a3, 0x2a9, 0x3a0, 0xd30,
    0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c, 0x53c, 0x435, 0x73f, 0x636,
    0x13a, 0x33,  0x339, 0x230, 0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895,
    0x99c, 0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99,  0x190, 0xf00, 0xe09,
    0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c, 0x70c, 0x605, 0x50f, 0x406, 0x30a,
    0x203, 0x109, 0x0};
static signed char triTable[256][16] = {
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
    {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
    {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
    {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
    {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
    {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
    {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
    {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
    {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
    {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
    {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
    {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
    {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
    {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
    {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
    {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
    {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
    {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
    {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
    {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
    {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
    {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
    {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
    {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
    {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
    {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
    {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
    {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
    {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
    {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
    {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
    {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
    {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
    {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
    {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
    {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
    {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
    {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
    {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
    {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
    {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
    {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
    {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
    {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
    {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
    {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
    {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
    {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
    {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
    {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
    {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
    {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
    {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
    {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
    {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
    {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
    {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
    {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
    {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
    {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
    {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
    {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
    {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
    {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
    {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
    {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
    {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
    {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
    {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
    {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
    {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
    {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
    {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
    {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
    {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
    {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
    {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
    {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
    {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
    {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
    {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
    {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
    {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
    {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
    {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
    {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
    {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
    {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
    {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
    {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
    {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
    {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
    {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
    {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
    {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
    {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
    {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
    {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
    {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
    {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
    {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
    {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
    {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
    {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
    {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
    {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
    {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
    {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
    {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
    {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
    {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
    {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
    {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
    {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
    {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
    {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
    {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
    {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
    {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
    {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
    {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
    {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
    {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
    {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
    {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
    {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
    {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
    {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
    {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
    {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
    {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
    {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
    {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
    {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
    {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
    {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
    {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
    {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
    {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
    {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
    {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
    {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
    {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
    {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
    {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
    {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
    {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
    {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
    {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
    {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};
//--------------------------------------------------------------------------------------------------------

class TriLinPoly {
    /* Compute a tri-linear polynomial within a given cube (i,j,k) x (i+1,j+1,k+1)
	 * Values are provided at the corners in CubeValues
	 * x,y,z must be defined on [0,1] where the length of the cube edge is one
	 */
    int ic, jc, kc;
    double a, b, c, d, e, f, g, h;
    double x, y, z;

public:
    DoubleArray Corners;

    TriLinPoly() { Corners.resize(2, 2, 2); }
    ~TriLinPoly() {}
    // Assign the polynomial within a cube from a mesh
    void assign(DoubleArray &A, int i, int j, int k) {
        // Save the cube indices
        ic = i;
        jc = j;
        kc = k;

        // local copy of the cube values
        Corners(0, 0, 0) = A(i, j, k);
        Corners(1, 0, 0) = A(i + 1, j, k);
        Corners(0, 1, 0) = A(i, j + 1, k);
        Corners(1, 1, 0) = A(i + 1, j + 1, k);
        Corners(0, 0, 1) = A(i, j, k + 1);
        Corners(1, 0, 1) = A(i + 1, j, k + 1);
        Corners(0, 1, 1) = A(i, j + 1, k + 1);
        Corners(1, 1, 1) = A(i + 1, j + 1, k + 1);

        // coefficients of the tri-linear approximation
        a = Corners(0, 0, 0);
        b = Corners(1, 0, 0) - a;
        c = Corners(0, 1, 0) - a;
        d = Corners(0, 0, 1) - a;
        e = Corners(1, 1, 0) - a - b - c;
        f = Corners(1, 0, 1) - a - b - d;
        g = Corners(0, 1, 1) - a - c - d;
        h = Corners(1, 1, 1) - a - b - c - d - e - f - g;
    }

    // Assign polynomial based on values manually assigned to Corners
    void assign() {
        ic = 0;
        jc = 0;
        kc = 0;

        // coefficients of the tri-linear approximation
        a = Corners(0, 0, 0);
        b = Corners(1, 0, 0) - a;
        c = Corners(0, 1, 0) - a;
        d = Corners(0, 0, 1) - a;
        e = Corners(1, 1, 0) - a - b - c;
        f = Corners(1, 0, 1) - a - b - d;
        g = Corners(0, 1, 1) - a - c - d;
        h = Corners(1, 1, 1) - a - b - c - d - e - f - g;
    }

    // Evaluate the polynomial at a point
    double eval(Point P) {
        double returnValue;
        x = P.x - 1.0 * ic;
        y = P.y - 1.0 * jc;
        z = P.z - 1.0 * kc;
        returnValue = a + b * x + c * y + d * z + e * x * y + f * x * z +
                      g * y * z + h * x * y * z;
        return returnValue;
    }
};
//--------------------------------------------------------------------------------------------------
/*
inline DoubleArray SOLVE( DoubleArray &A, DoubleArray &b)
{
	// solves the system A*x = b exactly

	DoubleArray AA = A.Copy();
	long int n = A.n;
	long int nrhs =1;
	IntArray piv(2*n);
	DoubleArray solution = b.Copy();

	long int info = 0;
	int ret = dgesv_(&n,&nrhs,AA.Pointer(),&n,(long int *)piv.Pointer(),solution.Pointer(),&n,&info);

	//	DTErrorMessage("Computation", string("Error=")+ DTInt2String(ret) + "info = " + DTInt2String(info));

    return solution;
}
//-----------------------------------------------------------------------------
inline Point NEWTON(Point x0,  DoubleArray &F, double &v,DoubleArray &S, int i, int j, int k, int &newton_steps)
{
	Point pt;
	// Performs a Newton iteration to compute the point x from initial guess x0

	double tol = 0.001;
	double dist = 1;

	// Map x0 into unit coordinates
	Point X0 = x0;
	x0.x = x0.x - i;
	x0.y = x0.y - j;
	x0.z = x0.z - k;

	// Compute coefficients for tri-linear approximations to the functions F and S
	// coefficients for F
	double Af, Bf, Cf, Df, Ef, Ff, Gf, Hf;
	Af = F(i,j,k);
	Bf = F(i+1,j,k)-Af;
	Cf = F(i,j+1,k)-Af;
	Df = F(i,j,k+1)-Af;
	Ef = F(i+1,j+1,k)-Af-Bf-Cf;
	Ff = F(i+1,j,k+1)-Af-Bf-Df;
	Gf = F(i,j+1,k+1)-Af-Cf-Df;
	Hf = F(i+1,j+1,k+1)-Af-Bf-Cf-Df-Ef-Ff-Gf;
	// coefficients for S
	double As, Bs, Cs, Ds, Es, Fs, Gs, Hs;
	As = S(i,j,k);
	Bs = S(i+1,j,k)-As;
	Cs = S(i,j+1,k)-As;
	Ds = S(i,j,k+1)-As;
	Es = S(i+1,j+1,k)-As-Bs-Cs;
	Fs = S(i+1,j,k+1)-As-Bs-Ds;
	Gs = S(i,j+1,k+1)-As-Cs-Ds;
	Hs = S(i+1,j+1,k+1)-As-Bs-Cs-Ds-Es-Fs-Gs;

	// Compute coefficients for plane in direction of grad F, grad S at point x0
	// Compute grad u = grad F , v = grad S
	int count; int number;
	count = 0;
	number = int(floor(newton_steps));
	//number = 2;
	pt = x0;

	bool pt_on_cube_leg = 0;

	// Compute the u = grad F, v = grad S and the normal vector N = u x v
	double u1, u2, u3, v1, v2, v3;
	u1 = Bf+Ef*x0.y+Ff*x0.z+Hf*x0.y*x0.z;
	u2 = Cf+Ef*x0.x+Gf*x0.z+Hf*x0.x*x0.z;
	u3 = Df+Ff*x0.x+Gf*x0.y+Hf*x0.x*x0.y;
	v1 = Bs+Es*x0.y+Fs*x0.z+Hs*x0.y*x0.z;
	v2 = Cs+Es*x0.x+Gs*x0.z+Hs*x0.x*x0.z;
	v3 = Ds+Fs*x0.x+Gs*x0.y+Hs*x0.x*x0.y;
	// Compute the normal vector N
	Point N;
	if (x0.x == 0 || x0.x == 1){
		if (x0.y == 0 || x0.y == 1){
			// on a cube leg
			pt_on_cube_leg = 1;
		}
		else if (x0.z == 0 || x0.z == 1){
			// on a cube leg
			pt_on_cube_leg = 1;
		}
		else {
			// use this cube face as the normal 
			N.x = 1;
			N.y = 0;
			N.z = 0;
		}
	}
	else if (x0.y == 0 || x0.y == 1){
		if (x0.z == 0 || x0.z == 1){
			// on a cube leg
			pt_on_cube_leg = 1;
		}
		else {
			// use this cube face as the normal
			N.x = 0;
			N.y = 1;
			N.z = 0;
		}
	}
	else if (x0.z == 0 || x0.z == 1){
		// use this cube face as the normal
		N.x = 0;
		N.y = 0;
		N.z = 1;
	}
	else {
		// Compute N = u x v
		N.x = u2*v3 - u3*v2;
		N.y = u3*v1 - u1*v3;
		N.z = u1*v2 - u2*v1;
	}

	if ( pt_on_cube_leg == 1){
		// return x0
		pt = x0;
	}
	else {
		while (count < number ){

			// Construct the Jacobean Matrix J
			DoubleArray J(3,3);
			J(0,0) = u1;
			J(0,1) = u2;
			J(0,2) = u3;
			J(1,0) = v1;
			J(1,1) = v2;
			J(1,2) = v3;
			J(2,0) = N.x;
			J(2,1) = N.y;
			J(2,2) = N.z;
			// Evaluate the vector M(x0)
			DoubleArray M(3);
			M(0) = -(Af+Bf*x0.x+Cf*x0.y+Df*x0.z+Ef*x0.x*x0.y+Ff*x0.x*x0.z+Gf*x0.y*x0.z+Hf*x0.x*x0.y*x0.z - v);
			M(1) = -(As+Bs*x0.x+Cs*x0.y+Ds*x0.z+Es*x0.x*x0.y+Fs*x0.x*x0.z+Gs*x0.y*x0.z+Hs*x0.x*x0.y*x0.z);
			M(2) = 0;
			// Compute the update to x0 using Newton's method
			Point dX;
			DoubleArray d = SOLVE(J,M);
			dX.x = d(0);
			dX.y = d(1);
			dX.z = d(2);

			if ( dX.x > 0.1 || dX.y > 0.1 || dX.z > 0.1 ||
				 dX.x < -0.1 || dX.y < -0.1 || dX.z < -0.1){
				count = number;
			}
			else {
				pt = x0+dX;
				x0 = pt;
			}
			M(0) = -(Af+Bf*x0.x+Cf*x0.y+Df*x0.z+Ef*x0.x*x0.y+Ff*x0.x*x0.z+Gf*x0.y*x0.z+Hf*x0.x*x0.y*x0.z - v);
			M(1) = -(As+Bs*x0.x+Cs*x0.y+Ds*x0.z+Es*x0.x*x0.y+Fs*x0.x*x0.z+Gs*x0.y*x0.z+Hs*x0.x*x0.y*x0.z);
			M(2) = 0;

			dist = sqrt(M(0)*M(0)+M(1)*M(1)+M(2)*M(2));

			// Compute u, v at new x0

			u1 = Bf+Ef*x0.y+Ff*x0.z+Hf*x0.y*x0.z;
			u2 = Cf+Ef*x0.x+Gf*x0.z+Hf*x0.x*x0.z;
			u3 = Df+Ff*x0.x+Gf*x0.y+Hf*x0.x*x0.y;
			v1 = Bs+Es*x0.y+Fs*x0.z+Hs*x0.y*x0.z;
			v2 = Cs+Es*x0.x+Gs*x0.z+Hs*x0.x*x0.z;
			v3 = Ds+Fs*x0.x+Gs*x0.y+Hs*x0.x*x0.y;

			count++;
		}
	}
	// map pt back to original coordinates
	pt.x = pt.x+i;
	pt.y = pt.y+j;
	pt.z = pt.z+k;

	return pt;
}
 */
//--------------------------------------------------------------------------------------------------

inline bool vertexcheck(Point &P, int n, int pos,
                        DTMutableList<Point> &cellvertices) {

    // returns true if P is a new vertex (one previously unencountered
    bool V = 1;
    int i = pos - n;
    for (i = pos - n; i < pos; i++) {
        if (P == cellvertices(i)) {
            V = 0;
        }
    }

    return V;
}

//--------------------------------------------------------------------------------------------------

inline bool ShareSide(Point &A, Point &B) {
    // returns true if points A and B share an x,y, or z coordinate
    bool l = 0;
    if (A.x != B.x && A.y != B.y && A.z != B.z) {
        l = 0;
    } else {
        if (floor(A.x) == A.x && floor(B.x) == B.x && A.x == B.x) {
            l = 1;
        }
        if (floor(A.y) == A.y && floor(B.y) == B.y && A.y == B.y) {
            l = 1;
        }
        if (floor(A.z) == A.z && floor(B.z) == B.z && A.z == B.z) {
            l = 1;
        }
    }

    return l;
}
//--------------------------------------------------------------------------------------------------

inline bool Interface(DoubleArray &A, const double v, int i, int j, int k) {
    // returns true if grid cell i, j, k contains a section of the interface
    bool Y = 0;

    if ((A(i, j, k) - v) * (A(i + 1, j, k) - v) < 0) {
        Y = 1;
    }
    // 2
    if ((A(i + 1, j, k) - v) * (A(i + 1, j + 1, k) - v) < 0) {
        Y = 1;
    }
    //3
    if ((A(i + 1, j + 1, k) - v) * (A(i, j + 1, k) - v) < 0) {
        Y = 1;
    }
    //4
    if ((A(i, j + 1, k) - v) * (A(i, j, k) - v) < 0) {
        Y = 1;
    }
    //5
    if ((A(i, j, k) - v) * (A(i, j, k + 1) - v) < 0) {
        Y = 1;
    }
    //6
    if ((A(i + 1, j, k) - v) * (A(i + 1, j, k + 1) - v) < 0) {
        Y = 1;
    }
    //7
    if ((A(i + 1, j + 1, k) - v) * (A(i + 1, j + 1, k + 1) - v) < 0) {
        Y = 1;
    }
    //8
    if ((A(i, j + 1, k) - v) * (A(i, j + 1, k + 1) - v) < 0) {
        Y = 1;
    }
    //9
    if ((A(i, j, k + 1) - v) * (A(i + 1, j, k + 1) - v) < 0) {
        Y = 1;
    }
    //10
    if ((A(i + 1, j, k + 1) - v) * (A(i + 1, j + 1, k + 1) - v) < 0) {
        Y = 1;
    }
    //11
    if ((A(i + 1, j + 1, k + 1) - v) * (A(i, j + 1, k + 1) - v) < 0) {
        Y = 1;
    }
    //12
    if ((A(i, j + 1, k + 1) - v) * (A(i, j, k + 1) - v) < 0) {
        Y = 1;
    }
    return Y;
}
//--------------------------------------------------------------------------------------------------

inline bool Fluid_Interface(DoubleArray &A, DoubleArray &S, const double v,
                            int i, int j, int k) {
    // returns true if grid cell i, j, k contains a section of the interface
    bool Y = 0;

    if ((A(i, j, k) - v) * (A(i + 1, j, k) - v) < 0 && S(i, j, k) > 0 &&
        S(i + 1, j, k) > 0) {
        Y = 1;
    }
    // 2
    if ((A(i + 1, j, k) - v) * (A(i + 1, j + 1, k) - v) < 0 &&
        S(i + 1, j, k) > 0 && S(i + 1, j + 1, k) > 0) {
        Y = 1;
    }
    //3
    if ((A(i + 1, j + 1, k) - v) * (A(i, j + 1, k) - v) < 0 &&
        S(i + 1, j + 1, k) > 0 && S(i, j + 1, k) > 0) {
        Y = 1;
    }
    //4
    if ((A(i, j + 1, k) - v) * (A(i, j, k) - v) < 0 && S(i, j, k) > 0 &&
        S(i, j + 1, k) > 0) {
        Y = 1;
    }
    //5
    if ((A(i, j, k) - v) * (A(i, j, k + 1) - v) < 0 && S(i, j, k) > 0 &&
        S(i, j, k + 1) > 0) {
        Y = 1;
    }
    //6
    if ((A(i + 1, j, k) - v) * (A(i + 1, j, k + 1) - v) < 0 &&
        S(i + 1, j, k) > 0 && S(i + 1, j, k + 1) > 0) {
        Y = 1;
    }
    //7
    if ((A(i + 1, j + 1, k) - v) * (A(i + 1, j + 1, k + 1) - v) < 0 &&
        S(i + 1, j + 1, k) > 0 && S(i + 1, j + 1, k + 1) > 0) {
        Y = 1;
    }
    //8
    if ((A(i, j + 1, k) - v) * (A(i, j + 1, k + 1) - v) < 0 &&
        S(i, j + 1, k) > 0 && S(i, j + 1, k + 1) > 0) {
        Y = 1;
    }
    //9
    if ((A(i, j, k + 1) - v) * (A(i + 1, j, k + 1) - v) < 0 &&
        S(i, j, k + 1) > 0 && S(i + 1, j, k + 1) > 0) {
        Y = 1;
    }
    //10
    if ((A(i + 1, j, k + 1) - v) * (A(i + 1, j + 1, k + 1) - v) < 0 &&
        S(i + 1, j, k + 1) > 0 && S(i + 1, j + 1, k + 1) > 0) {
        Y = 1;
    }
    //11
    if ((A(i + 1, j + 1, k + 1) - v) * (A(i, j + 1, k + 1) - v) < 0 &&
        S(i + 1, j + 1, k + 1) > 0 && S(i, j + 1, k + 1) > 0) {
        Y = 1;
    }
    //12
    if ((A(i, j + 1, k + 1) - v) * (A(i, j, k + 1) - v) < 0 &&
        S(i, j, k + 1) > 0 && S(i, j + 1, k + 1) > 0) {
        Y = 1;
    }
    return Y;
}
//--------------------------------------------------------------------------------------------------
inline bool Solid(DoubleArray &A, int i, int j, int k) {

    bool X = 0;

    // return 0 if there is no solid phase in grid cell i,j,k

    if (A(i, j, k) == 0) {
        X = 1;
    }
    if (A(i + 1, j, k) == 0) {
        X = 1;
    }
    if (A(i, j + 1, k) == 0) {
        X = 1;
    }
    if (A(i, j, k + 1) == 0) {
        X = 1;
    }
    if (A(i + 1, j + 1, k) == 0) {
        X = 1;
    }
    if (A(i + 1, j, k + 1) == 0) {
        X = 1;
    }
    if (A(i, j + 1, k + 1) == 0) {
        X = 1;
    }
    if (A(i + 1, j + 1, k + 1) == 0) {
        X = 1;
    }
    return X;
}
//-------------------------------------------------------------------------------
inline Point VertexInterp(const Point &p1, const Point &p2, double valp1,
                          double valp2) {
    return (p1 + (-valp1 / (valp2 - valp1)) * (p2 - p1));
}
//-------------------------------------------------------------------------------
inline void SolidMarchingCubes(DoubleArray &A, const double &v, DoubleArray &B,
                               const double &isovalue, int i, int j, int k,
                               int m, int n, int o,
                               DTMutableList<Point> &cellvertices,
                               int &lengthvertices, IntArray &Triangles,
                               int &nTris, DoubleArray &values) {

    NULL_USE(isovalue);
    NULL_USE(m);
    NULL_USE(n);
    NULL_USE(o);

    // THIS SUBROUTINE COMPUTES THE VERTICES FOR THE SOLID PHASE IN
    // A PARTICULAR GRID CELL, THEN ARRANGES THEM INTO TRIANGLES
    // ALSO ORGANIZES THE LIST OF VALUES TO CORRESPOND WITH VERTEX LIST

    NULL_USE(v);

    //int N = 0;
    Point P, Q;
    Point PlaceHolder;
    //int pos = lengthvertices;
    double temp;
    Point C0, C1, C2, C3, C4, C5, C6, C7;

    int TriangleCount;
    int NewVertexCount;
    int CubeIndex;

    Point VertexList[12];
    Point NewVertexList[12];
    double ValueList[12];
    double NewValueList[12];
    int LocalRemap[12];

    // int m; int n; int o;
    //int p; int q;

    // Values from array 'A' at the cube corners
    double CubeValues[8];
    CubeValues[0] = A(i, j, k);
    CubeValues[1] = A(i + 1, j, k);
    CubeValues[2] = A(i + 1, j + 1, k);
    CubeValues[3] = A(i, j + 1, k);
    CubeValues[4] = A(i, j, k + 1);
    CubeValues[5] = A(i + 1, j, k + 1);
    CubeValues[6] = A(i + 1, j + 1, k + 1);
    CubeValues[7] = A(i, j + 1, k + 1);

    // Values from array 'B' at the cube corners
    double CubeValues_B[8];
    CubeValues_B[0] = B(i, j, k);
    CubeValues_B[1] = B(i + 1, j, k);
    CubeValues_B[2] = B(i + 1, j + 1, k);
    CubeValues_B[3] = B(i, j + 1, k);
    CubeValues_B[4] = B(i, j, k + 1);
    CubeValues_B[5] = B(i + 1, j, k + 1);
    CubeValues_B[6] = B(i + 1, j + 1, k + 1);
    CubeValues_B[7] = B(i, j + 1, k + 1);

    // Points corresponding to cube corners
    C0.x = 0.0;
    C0.y = 0.0;
    C0.z = 0.0;
    C1.x = 1.0;
    C1.y = 0.0;
    C1.z = 0.0;
    C2.x = 1.0;
    C2.y = 1.0;
    C2.z = 0.0;
    C3.x = 0.0;
    C3.y = 1.0;
    C3.z = 0.0;
    C4.x = 0.0;
    C4.y = 0.0;
    C4.z = 1.0;
    C5.x = 1.0;
    C5.y = 0.0;
    C5.z = 1.0;
    C6.x = 1.0;
    C6.y = 1.0;
    C6.z = 1.0;
    C7.x = 0.0;
    C7.y = 1.0;
    C7.z = 1.0;
    //Determine the index into the edge table which
    //tells us which vertices are inside of the surface
    CubeIndex = 0;
    if (CubeValues[0] < 0.0f)
        CubeIndex |= 1;
    if (CubeValues[1] < 0.0f)
        CubeIndex |= 2;
    if (CubeValues[2] < 0.0f)
        CubeIndex |= 4;
    if (CubeValues[3] < 0.0f)
        CubeIndex |= 8;
    if (CubeValues[4] < 0.0f)
        CubeIndex |= 16;
    if (CubeValues[5] < 0.0f)
        CubeIndex |= 32;
    if (CubeValues[6] < 0.0f)
        CubeIndex |= 64;
    if (CubeValues[7] < 0.0f)
        CubeIndex |= 128;

    //Find the vertices where the surface intersects the cube
    if (edgeTable[CubeIndex] & 1) {
        P = VertexInterp(C0, C1, CubeValues[0], CubeValues[1]);
        VertexList[0] = P;
        Q = C0;
        temp = sqrt((P.x - Q.x) * (P.x - Q.x) + (P.y - Q.y) * (P.y - Q.y) +
                    (P.z - Q.z) * (P.z - Q.z));
        ValueList[0] =
            CubeValues_B[0] + temp * (CubeValues_B[1] - CubeValues_B[0]);
    }
    if (edgeTable[CubeIndex] & 2) {
        P = VertexInterp(C1, C2, CubeValues[1], CubeValues[2]);
        VertexList[1] = P;
        Q = C1;
        temp = sqrt((P.x - Q.x) * (P.x - Q.x) + (P.y - Q.y) * (P.y - Q.y) +
                    (P.z - Q.z) * (P.z - Q.z));
        ValueList[1] =
            CubeValues_B[1] + temp * (CubeValues_B[2] - CubeValues_B[1]);
    }
    if (edgeTable[CubeIndex] & 4) {
        P = VertexInterp(C2, C3, CubeValues[2], CubeValues[3]);
        VertexList[2] = P;
        Q = C2;
        temp = sqrt((P.x - Q.x) * (P.x - Q.x) + (P.y - Q.y) * (P.y - Q.y) +
                    (P.z - Q.z) * (P.z - Q.z));
        ValueList[2] =
            CubeValues_B[2] + temp * (CubeValues_B[3] - CubeValues_B[2]);
    }
    if (edgeTable[CubeIndex] & 8) {
        P = VertexInterp(C3, C0, CubeValues[3], CubeValues[0]);
        VertexList[3] = P;
        Q = C3;
        temp = sqrt((P.x - Q.x) * (P.x - Q.x) + (P.y - Q.y) * (P.y - Q.y) +
                    (P.z - Q.z) * (P.z - Q.z));
        ValueList[3] =
            CubeValues_B[3] + temp * (CubeValues_B[0] - CubeValues_B[3]);
    }
    if (edgeTable[CubeIndex] & 16) {
        P = VertexInterp(C4, C5, CubeValues[4], CubeValues[5]);
        VertexList[4] = P;
        Q = C4;
        temp = sqrt((P.x - Q.x) * (P.x - Q.x) + (P.y - Q.y) * (P.y - Q.y) +
                    (P.z - Q.z) * (P.z - Q.z));
        ValueList[4] =
            CubeValues_B[4] + temp * (CubeValues_B[5] - CubeValues_B[4]);
    }
    if (edgeTable[CubeIndex] & 32) {
        P = VertexInterp(C5, C6, CubeValues[5], CubeValues[6]);
        VertexList[5] = P;
        Q = C5;
        temp = sqrt((P.x - Q.x) * (P.x - Q.x) + (P.y - Q.y) * (P.y - Q.y) +
                    (P.z - Q.z) * (P.z - Q.z));
        ValueList[5] =
            CubeValues_B[5] + temp * (CubeValues_B[6] - CubeValues_B[5]);
    }
    if (edgeTable[CubeIndex] & 64) {
        P = VertexInterp(C6, C7, CubeValues[6], CubeValues[7]);
        VertexList[6] = P;
        Q = C6;
        temp = sqrt((P.x - Q.x) * (P.x - Q.x) + (P.y - Q.y) * (P.y - Q.y) +
                    (P.z - Q.z) * (P.z - Q.z));
        ValueList[6] =
            CubeValues_B[6] + temp * (CubeValues_B[7] - CubeValues_B[6]);
    }
    if (edgeTable[CubeIndex] & 128) {
        P = VertexInterp(C7, C4, CubeValues[7], CubeValues[4]);
        VertexList[7] = P;
        Q = C7;
        temp = sqrt((P.x - Q.x) * (P.x - Q.x) + (P.y - Q.y) * (P.y - Q.y) +
                    (P.z - Q.z) * (P.z - Q.z));
        ValueList[7] =
            CubeValues_B[7] + temp * (CubeValues_B[4] - CubeValues_B[7]);
    }
    if (edgeTable[CubeIndex] & 256) {
        P = VertexInterp(C0, C4, CubeValues[0], CubeValues[4]);
        VertexList[8] = P;
        Q = C0;
        temp = sqrt((P.x - Q.x) * (P.x - Q.x) + (P.y - Q.y) * (P.y - Q.y) +
                    (P.z - Q.z) * (P.z - Q.z));
        ValueList[8] =
            CubeValues_B[0] + temp * (CubeValues_B[4] - CubeValues_B[0]);
    }
    if (edgeTable[CubeIndex] & 512) {
        P = VertexInterp(C1, C5, CubeValues[1], CubeValues[5]);
        VertexList[9] = P;
        Q = C1;
        temp = sqrt((P.x - Q.x) * (P.x - Q.x) + (P.y - Q.y) * (P.y - Q.y) +
                    (P.z - Q.z) * (P.z - Q.z));
        ValueList[9] =
            CubeValues_B[1] + temp * (CubeValues_B[5] - CubeValues_B[1]);
    }
    if (edgeTable[CubeIndex] & 1024) {
        P = VertexInterp(C2, C6, CubeValues[2], CubeValues[6]);
        VertexList[10] = P;
        Q = C2;
        temp = sqrt((P.x - Q.x) * (P.x - Q.x) + (P.y - Q.y) * (P.y - Q.y) +
                    (P.z - Q.z) * (P.z - Q.z));
        ValueList[10] =
            CubeValues_B[2] + temp * (CubeValues_B[6] - CubeValues_B[2]);
    }
    if (edgeTable[CubeIndex] & 2048) {
        P = VertexInterp(C3, C7, CubeValues[3], CubeValues[7]);
        VertexList[11] = P;
        Q = C3;
        temp = sqrt((P.x - Q.x) * (P.x - Q.x) + (P.y - Q.y) * (P.y - Q.y) +
                    (P.z - Q.z) * (P.z - Q.z));
        ValueList[11] =
            CubeValues_B[3] + temp * (CubeValues_B[7] - CubeValues_B[3]);
    }

    NewVertexCount = 0;
    for (int idx = 0; idx < 12; idx++)
        LocalRemap[idx] = -1;

    for (int idx = 0; triTable[CubeIndex][idx] != -1; idx++) {
        if (LocalRemap[triTable[CubeIndex][idx]] == -1) {
            NewVertexList[NewVertexCount] =
                VertexList[triTable[CubeIndex][idx]];
            NewValueList[NewVertexCount] = ValueList[triTable[CubeIndex][idx]];
            LocalRemap[triTable[CubeIndex][idx]] = NewVertexCount;
            NewVertexCount++;
        }
    }

    for (int idx = 0; idx < NewVertexCount; idx++) {
        P = NewVertexList[idx];
        P.x += i;
        P.y += j;
        P.z += k;
        cellvertices(idx) = P;
        values(idx) = NewValueList[idx];
    }
    lengthvertices = NewVertexCount;

    TriangleCount = 0;
    for (int idx = 0; triTable[CubeIndex][idx] != -1; idx += 3) {
        Triangles(0, TriangleCount) = LocalRemap[triTable[CubeIndex][idx + 0]];
        Triangles(1, TriangleCount) = LocalRemap[triTable[CubeIndex][idx + 1]];
        Triangles(2, TriangleCount) = LocalRemap[triTable[CubeIndex][idx + 2]];
        TriangleCount++;
    }
    nTris = TriangleCount;

    // *    *    *   ARRANGE VERTICES SO THAT NEIGHBORS SHARE A FACE    *    *    *
    // *    *    *    PERFORM SAME OPERATIONS TO THE LIST OF VALUES     *    *    *

    /*	pos = N = lengthvertices;
	
    for (q = pos-N; q < pos-2; q++) {
        for (o = q+2; o < pos-1; o++) {
            if (ShareSide(cellvertices(q), cellvertices(o)) == 1) {
                PlaceHolder = cellvertices(q+1);
                cellvertices(q+1) = cellvertices(o);
                cellvertices(o) = PlaceHolder;
				
				temp = values(q+1);
				values(q+1) = values(o);
				values(o) = temp;
				// Swap the triangle lists as needed
				for (int idx=0; idx<TriangleCount; idx++){
					if (Triangles(0,idx) == q+1)	Triangles(0,idx) = o;
					else if (Triangles(0,idx) == o)	Triangles(0,idx) = q+1;
					if (Triangles(1,idx) == q+1)	Triangles(1,idx) = o;
					else if (Triangles(1,idx) == o)	Triangles(1,idx) = q+1;
					if (Triangles(2,idx) == q+1)	Triangles(2,idx) = o;
					else if (Triangles(2,idx) == o)	Triangles(2,idx) = q+1;
				}
            }
        }
		
		// make sure other neighbor of vertex 1 is in last spot
        if (q == pos-N){
            for (p = q+2; p < pos-1; p++){
                if (ShareSide(cellvertices(q), cellvertices(p)) == 1){
                    PlaceHolder = cellvertices(pos-1);
                    cellvertices(pos-1) = cellvertices(p);
                    cellvertices(p) = PlaceHolder;
					
					temp = values(pos-1);
					values(pos-1) = values(p);
					values(p) = temp;
					// Swap the triangle lists as needed
					for (int idx=0; idx<TriangleCount; idx++){
						if (Triangles(0,idx) == p)	Triangles(0,idx) = pos-1;
						else if (Triangles(0,idx) == pos-1)	Triangles(0,idx) = p;
						if (Triangles(1,idx) == p)	Triangles(1,idx) = pos-1;
						else if (Triangles(1,idx) == pos-1)	Triangles(1,idx) = p;
						if (Triangles(2,idx) == p)	Triangles(2,idx) = pos-1;
						else if (Triangles(2,idx) == pos-1)	Triangles(2,idx) = p;
					}

                }
            }
        }
        if ( ShareSide(cellvertices(pos-2), cellvertices(pos-3)) != 1 ){
            if (ShareSide( cellvertices(pos-3), cellvertices(pos-1)) == 1 && ShareSide( cellvertices(pos-N),
																					   cellvertices(pos-2)) == 1 ){
                PlaceHolder = cellvertices(pos-2);
                cellvertices(pos-2) = cellvertices(pos-1);
                cellvertices(pos-1) = PlaceHolder;
				
				temp = values(pos-2);
				values(pos-2) = values(pos-1);
				values(pos-1) = temp;
				
				// Swap the triangle lists as needed
				for (int idx=0; idx<TriangleCount; idx++){
					if (Triangles(0,idx) == pos-1)	Triangles(0,idx) = pos-2;
					else if (Triangles(0,idx) == pos-2)	Triangles(0,idx) = pos-1;
					if (Triangles(1,idx) == pos-1)	Triangles(1,idx) = pos-2;
					else if (Triangles(1,idx) == pos-2)	Triangles(1,idx) = pos-1;
					if (Triangles(2,idx) == pos-1)	Triangles(2,idx) = pos-2;
					else if (Triangles(2,idx) == pos-2)	Triangles(2,idx) = pos-1;
				}

            }
        }
        if ( ShareSide(cellvertices(pos-1), cellvertices(pos-2)) != 1 ){
            if (ShareSide( cellvertices(pos-3), cellvertices(pos-1)) == 1 && ShareSide(cellvertices(pos-4),
																					   cellvertices(pos-2)) == 1 ){
                PlaceHolder = cellvertices(pos-3);
                cellvertices(pos-3) = cellvertices(pos-2);
                cellvertices(pos-2) = PlaceHolder;
				
				temp = values(pos-3);
				values(pos-3) = values(pos-2);
				values(pos-2) = temp;
				
				// Swap the triangle lists as needed
				for (int idx=0; idx<TriangleCount; idx++){
					if (Triangles(0,idx) == pos-2)	Triangles(0,idx) = pos-3;
					else if (Triangles(0,idx) == pos-2)	Triangles(0,idx) = pos-2;
					if (Triangles(1,idx) == pos-2)	Triangles(1,idx) = pos-3;
					else if (Triangles(1,idx) == pos-2)	Triangles(1,idx) = pos-2;
					if (Triangles(2,idx) == pos-2)	Triangles(2,idx) = pos-3;
					else if (Triangles(2,idx) == pos-2)	Triangles(2,idx) = pos-2;
				}
            }
            if (ShareSide( cellvertices(pos-N+1), cellvertices(pos-3)) == 1 && ShareSide(cellvertices(pos-1),
																						 cellvertices(pos-N+1)) == 1 ){
                PlaceHolder = cellvertices(pos-2);
                cellvertices(pos-2) = cellvertices(pos-N+1);
                cellvertices(pos-N+1) = PlaceHolder;
				
				temp = values(pos-2);
				values(pos-2) = values(pos-N+1);
				values(pos-N+1) = temp;
				
				// Swap the triangle lists as needed
				for (int idx=0; idx<TriangleCount; idx++){
					if (Triangles(0,idx) == pos-N+1)	Triangles(0,idx) = pos-2;
					else if (Triangles(0,idx) == pos-2)	Triangles(0,idx) = pos-N+1;
					if (Triangles(1,idx) == pos-N+1)	Triangles(1,idx) = pos-2;
					else if (Triangles(1,idx) == pos-2)	Triangles(1,idx) = pos-N+1;
					if (Triangles(2,idx) == pos-N+1)	Triangles(2,idx) = pos-2;
					else if (Triangles(2,idx) == pos-2)	Triangles(2,idx) = pos-N+1;
				}
            }
        }
        if ( ShareSide(cellvertices(pos-N), cellvertices(pos-N+1)) != 1 ){
            if (ShareSide( cellvertices(pos-N), cellvertices(pos-2)) == 1 && ShareSide(cellvertices(pos-1),
																					   cellvertices(pos-N+1)) == 1){
                PlaceHolder = cellvertices(pos-1);
                cellvertices(pos-1) = cellvertices(pos-N);
                cellvertices(pos-N) = PlaceHolder;
				
				temp = values(pos-1);
				values(pos-1) = values(pos-N);
				values(pos-N) = temp;
				
				// Swap the triangle lists as needed
				for (int idx=0; idx<TriangleCount; idx++){
					if (Triangles(0,idx) == pos-N)	Triangles(0,idx) = pos-1;
					else if (Triangles(0,idx) == pos-1)	Triangles(0,idx) = pos-N;
					if (Triangles(1,idx) == pos-N)	Triangles(1,idx) = pos-1;
					else if (Triangles(1,idx) == pos-1)	Triangles(1,idx) = pos-N;
					if (Triangles(2,idx) == pos-N)	Triangles(2,idx) = pos-1;
					else if (Triangles(2,idx) == pos-1)	Triangles(2,idx) = pos-N;
				}
            }
        }
    }
 */
}
//-------------------------------------------------------------------------------
inline void SOL_SURF(DoubleArray &A, const double &v, DoubleArray &B,
                     const double &isovalue, int i, int j, int k, int m, int n,
                     int o, DTMutableList<Point> &cellvertices,
                     int &lengthvertices, IntArray &Tlist, int &nTris,
                     DoubleArray &values) {
    NULL_USE(isovalue);
    NULL_USE(m);
    NULL_USE(n);
    NULL_USE(o);

    // THIS SUBROUTINE COMPUTES THE VERTICES FOR THE SOLID PHASE IN
    // A PARTICULAR GRID CELL, THEN ARRANGES THEM INTO TRIANGLES
    // ALSO ORGANIZES THE LIST OF VALUES TO CORRESPOND WITH VERTEX LIST

    int N = 0;
    Point P;
    Point PlaceHolder;
    int pos = lengthvertices;
    float temp;

    // int m; int n; int o;
    int p;
    int q;

    // for each leg of the triangle, see if a vertex exists
    // if so, find the vertex, and perform an extrapolation

    // Go over each corner -- check to see if the corners are themselves vertices
    //1
    if (A(i, j, k) == v) {
        P.x = i;
        P.y = j;
        P.z = k;
        values(pos) = B(i, j, k);
        cellvertices(pos++) = P;
        N++;
    }
    //2
    if (A(i + 1, j, k) == v) {
        P.x = i + 1;
        P.y = j;
        P.z = k;
        values(pos) = B(i + 1, j, k);
        cellvertices(pos++) = P;
        N++;
    }
    //3
    if (A(i + 1, j + 1, k) == v) {
        P.x = i + 1;
        P.y = j + 1;
        P.z = k;
        values(pos) = B(i + 1, j + 1, k);
        cellvertices(pos++) = P;
        N++;
    }
    //4
    if (A(i, j + 1, k) == v) {
        P.x = i;
        P.y = j + 1;
        P.z = k;
        values(pos) = B(i, j + 1, k);
        cellvertices(pos++) = P;
        N++;
    }
    //5
    if (A(i, j, k + 1) == v) {
        P.x = i;
        P.y = j;
        P.z = k + 1;
        values(pos) = B(i, j, k + 1);
        cellvertices(pos++) = P;
        N++;
    }
    //6
    if (A(i + 1, j, k + 1) == v) {
        P.x = i + 1;
        P.y = j;
        P.z = k + 1;
        values(pos) = B(i + 1, j, k + 1);
        cellvertices(pos++) = P;
        N++;
    }
    //7
    if (A(i + 1, j + 1, k + 1) == v) {
        P.x = i + 1;
        P.y = j + 1;
        P.z = k + 1;
        values(pos) = B(i + 1, j + 1, k + 1);
        cellvertices(pos++) = P;
        N++;
    }
    //8
    if (A(i, j + 1, k + 1) == v) {
        P.x = i;
        P.y = j + 1;
        P.z = k + 1;
        values(pos) = B(i, j + 1, k + 1);
        cellvertices(pos++) = P;
        N++;
    }

    // Go through each side, compute P for sides of box spiraling up

    if ((A(i, j, k) - v) * (A(i + 1, j, k) - v) < 0) {
        P.x = i + (A(i, j, k) - v) / (A(i, j, k) - A(i + 1, j, k));
        P.y = j;
        P.z = k;
        // compute extrapolated fluid value at P
        //		if ( A(i,j,k) > v){
        //			 values(pos) = EXTRAP(B, isovalue, i,j,k,4, P);
        //		}
        //		else{
        //			values(pos) = EXTRAP(B,isovalue, i+1,j,k, 1, P);
        //		}
        // Interpolate value at P using padded mesh B
        values(pos) = B(i, j, k) * (1 - P.x + i) + B(i + 1, j, k) * (P.x - i);

        cellvertices(pos++) = P;
        N++;
    }
    // 2
    if ((A(i + 1, j, k) - v) * (A(i + 1, j + 1, k) - v) < 0) {
        P.x = i + 1;
        P.y = j + (A(i + 1, j, k) - v) / (A(i + 1, j, k) - A(i + 1, j + 1, k));
        P.z = k;
        if (vertexcheck(P, N, pos, cellvertices) ==
            1) { // P is a new vertex (not counted twice)
                 // compute extrapolated fluid value at P
                 //			if ( A(i+1,j,k) > v){
                 //				values(pos) = EXTRAP(B,isovalue, i+1,j,k, 5, P);
                 //			}
                 //			else{
                 //				values(pos) = EXTRAP(B,isovalue, i+1,j+1,k, 2, P);
                 //			}
                 // Interpolate value at P using padded mesh B
            values(pos) =
                B(i + 1, j, k) * (1 - P.y + j) + B(i + 1, j + 1, k) * (P.y - j);
            cellvertices(pos++) = P;
            N++;
        }
    }
    //3
    if ((A(i + 1, j + 1, k) - v) * (A(i, j + 1, k) - v) < 0) {
        P.x = i + (A(i, j + 1, k) - v) / (A(i, j + 1, k) - A(i + 1, j + 1, k));
        P.y = j + 1;
        P.z = k;
        if (vertexcheck(P, N, pos, cellvertices) ==
            1) { // P is a new vertex (not counted twice)
                 // compute extrapolated fluid value at P
                 //			if ( A(i+1,j+1,k) > v){
                 //				values(pos) = EXTRAP(B,isovalue, i+1,j+1,k, 1, P);
                 //			}
                 //			else{
                 //				values(pos) = EXTRAP(B,isovalue, i,j+1,k, 4, P);
                 //			}
                 // Interpolate value at P using padded mesh B
            values(pos) =
                B(i, j + 1, k) * (1 - P.x + i) + B(i + 1, j + 1, k) * (P.x - i);
            cellvertices(pos++) = P;
            N++;
        }
    }
    //4
    if ((A(i, j + 1, k) - v) * (A(i, j, k) - v) < 0) {
        P.x = i;
        P.y = j + (A(i, j, k) - v) / (A(i, j, k) - A(i, j + 1, k));
        P.z = k;
        if (vertexcheck(P, N, pos, cellvertices) ==
            1) { // P is a new vertex (not counted twice)
                 // compute extrapolated fluid value at P
                 //			if ( A(i,j,k) > v){
                 //				values(pos) = EXTRAP(B,isovalue, i,j,k, 5, P);
                 //			}
                 //			else{
                 //				values(pos) = EXTRAP(B,isovalue, i,j+1,k, 2, P);
                 //			}
                 // Interpolate value at P using padded mesh B
            values(pos) =
                B(i, j, k) * (1 - P.y + j) + B(i, j + 1, k) * (P.y - j);
            cellvertices(pos++) = P;
            N++;
        }
    }
    //5
    if ((A(i, j, k) - v) * (A(i, j, k + 1) - v) < 0) {
        P.x = i;
        P.y = j;
        P.z = k + (A(i, j, k) - v) / (A(i, j, k) - A(i, j, k + 1));
        if (vertexcheck(P, N, pos, cellvertices) ==
            1) { // P is a new vertex (not counted twice)
                 // compute extrapolated fluid value at P
                 //			if ( A(i,j,k) > v){
                 //				values(pos) = EXTRAP(B,isovalue, i,j,k, 6, P);
                 //			}
                 //			else{
                 //				values(pos) = EXTRAP(B,isovalue, i,j,k+1, 3, P);
                 //			}

            // Interpolate value at P using padded mesh B
            values(pos) =
                B(i, j, k) * (1 - P.z + k) + B(i, j, k + 1) * (P.z - k);
            cellvertices(pos++) = P;
            N++;
        }
    }
    //6
    if ((A(i + 1, j, k) - v) * (A(i + 1, j, k + 1) - v) < 0) {
        P.x = i + 1;
        P.y = j;
        P.z = k + (A(i + 1, j, k) - v) / (A(i + 1, j, k) - A(i + 1, j, k + 1));
        if (vertexcheck(P, N, pos, cellvertices) ==
            1) { // P is a new vertex (not counted twice)
                 // compute extrapolated fluid value at P
                 //			if ( A(i+1,j,k) > v){
                 //				values(pos) = EXTRAP(B,isovalue, i+1,j,k, 6, P);
                 //			}
                 //			else{
                 //				values(pos) = EXTRAP(B,isovalue, i+1,j,k+1, 3, P);
                 //			}

            // Interpolate value at P using padded mesh B
            values(pos) =
                B(i + 1, j, k) * (1 - P.z + k) + B(i + 1, j, k + 1) * (P.z - k);
            cellvertices(pos++) = P;
            N++;
        }
    }
    //7
    if ((A(i + 1, j + 1, k) - v) * (A(i + 1, j + 1, k + 1) - v) < 0) {
        P.x = i + 1;
        P.y = j + 1;
        P.z = k + (A(i + 1, j + 1, k) - v) /
                      (A(i + 1, j + 1, k) - A(i + 1, j + 1, k + 1));
        if (vertexcheck(P, N, pos, cellvertices) ==
            1) { // P is a new vertex (not counted twice)
                 // compute extrapolated fluid value at P
                 //			if ( A(i+1,j+1,k) > v){
                 //				values(pos) = EXTRAP(B,isovalue, i+1,j+1,k, 6, P);
                 //			}
                 //			else{
                 //				values(pos) = EXTRAP(B,isovalue, i+1,j+1,k+1, 3, P);
                 //			}

            // Interpolate value at P using padded mesh B
            values(pos) = B(i + 1, j + 1, k) * (1 - P.z + k) +
                          B(i + 1, j + 1, k + 1) * (P.z - k);
            cellvertices(pos++) = P;
            N++;
        }
    }
    //8
    if ((A(i, j + 1, k) - v) * (A(i, j + 1, k + 1) - v) < 0) {
        P.x = i;
        P.y = j + 1;
        P.z = k + (A(i, j + 1, k) - v) / (A(i, j + 1, k) - A(i, j + 1, k + 1));
        if (vertexcheck(P, N, pos, cellvertices) ==
            1) { // P is a new vertex (not counted twice)
                 // compute extrapolated fluid value at P
                 //			if ( A(i,j+1,k) > v){
                 //				values(pos) = EXTRAP(B,isovalue, i,j+1,k, 6, P);
                 //			}
                 //			else{
                 //				values(pos) = EXTRAP(B,isovalue, i,j+1,k+1, 3, P);
                 //			}

            // Interpolate value at P using padded mesh B
            values(pos) =
                B(i, j + 1, k) * (1 - P.z + k) + B(i, j + 1, k + 1) * (P.z - k);
            cellvertices(pos++) = P;
            N++;
        }
    }
    //9
    if ((A(i, j, k + 1) - v) * (A(i + 1, j, k + 1) - v) < 0) {
        P.x = i + (A(i, j, k + 1) - v) / (A(i, j, k + 1) - A(i + 1, j, k + 1));
        P.y = j;
        P.z = k + 1;
        if (vertexcheck(P, N, pos, cellvertices) ==
            1) { // P is a new vertex (not counted twice)
                 // compute extrapolated fluid value at P
                 //			if ( A(i,j,k+1) > v){
                 //				values(pos) = EXTRAP(B,isovalue, i,j,k+1, 4, P);
                 //			}
                 //			else{
                 //				values(pos) = EXTRAP(B,isovalue, i+1,j,k+1, 1, P);
                 //			}

            // Interpolate value at P using padded mesh B
            values(pos) =
                B(i, j, k + 1) * (1 - P.x + i) + B(i + 1, j, k + 1) * (P.x - i);
            cellvertices(pos++) = P;
            N++;
        }
    }
    //10
    if ((A(i + 1, j, k + 1) - v) * (A(i + 1, j + 1, k + 1) - v) < 0) {
        P.x = i + 1;
        P.y = j + (A(i + 1, j, k + 1) - v) /
                      (A(i + 1, j, k + 1) - A(i + 1, j + 1, k + 1));
        P.z = k + 1;
        if (vertexcheck(P, N, pos, cellvertices) ==
            1) { // P is a new vertex (not counted twice)
                 // compute extrapolated fluid value at P
                 //			if ( A(i+1,j,k+1) > v){
                 //				values(pos) = EXTRAP(B,isovalue, i+1,j,k+1, 5, P);
                 //			}
                 //			else{
                 //				values(pos) = EXTRAP(B,isovalue, i+1,j+1,k+1, 2, P);
                 //			}

            // Interpolate value at P using padded mesh B
            values(pos) = B(i + 1, j, k + 1) * (1 - P.y + j) +
                          B(i + 1, j + 1, k + 1) * (P.y - j);
            cellvertices(pos++) = P;
            N++;
        }
    }
    //11
    if ((A(i + 1, j + 1, k + 1) - v) * (A(i, j + 1, k + 1) - v) < 0) {
        P.x = i + (A(i, j + 1, k + 1) - v) /
                      (A(i, j + 1, k + 1) - A(i + 1, j + 1, k + 1));
        P.y = j + 1;
        P.z = k + 1;
        if (vertexcheck(P, N, pos, cellvertices) ==
            1) { // P is a new vertex (not counted twice)
                 // compute extrapolated fluid value at P
                 //			if ( A(i,j+1,k+1) > v){
                 //				values(pos) = EXTRAP(B,isovalue, i,j+1,k+1, 4, P);
                 //			}
                 //			else{
                 //				values(pos) = EXTRAP(B,isovalue, i+1,j+1,k+1, 1, P);
                 //			}

            // Interpolate value at P using padded mesh B
            values(pos) = B(i, j + 1, k + 1) * (1 - P.x + i) +
                          B(i + 1, j + 1, k + 1) * (P.x - i);
            cellvertices(pos++) = P;
            N++;
        }
    }
    //12
    if ((A(i, j + 1, k + 1) - v) * (A(i, j, k + 1) - v) < 0) {
        P.x = i;
        P.y = j + (A(i, j, k + 1) - v) / (A(i, j, k + 1) - A(i, j + 1, k + 1));
        P.z = k + 1;
        if (vertexcheck(P, N, pos, cellvertices) ==
            1) { // P is a new vertex (not counted twice)
                 // compute extrapolated fluid value at P
                 //			if ( A(i,j,k+1) > v){
                 //				values(pos) = EXTRAP(B,isovalue, i,j,k+1, 5, P);
                 //			}
                 //			else{
                 //				values(pos) = EXTRAP(B,isovalue, i,j+1,k+1, 2, P);
                 //			}

            // Interpolate value at P using padded mesh B
            values(pos) =
                B(i, j, k + 1) * (1 - P.y + j) + B(i, j + 1, k + 1) * (P.y - j);
            cellvertices(pos++) = P;
            N++;
        }
    }

    lengthvertices = pos;

    // *    *    *   ARRANGE VERTICES SO THAT NEIGHBORS SHARE A FACE    *    *    *
    // *    *    *    PERFORM SAME OPERATIONS TO THE LIST OF VALUES     *    *    *

    for (q = pos - N; q < pos - 2; q++) {
        for (o = q + 2; o < pos - 1; o++) {
            if (ShareSide(cellvertices(q), cellvertices(o)) == 1) {
                PlaceHolder = cellvertices(q + 1);
                cellvertices(q + 1) = cellvertices(o);
                cellvertices(o) = PlaceHolder;

                temp = values(q + 1);
                values(q + 1) = values(o);
                values(o) = temp;
            }
        }

        // make sure other neighbor of vertex 1 is in last spot
        if (q == pos - N) {
            for (p = q + 2; p < pos - 1; p++) {
                if (ShareSide(cellvertices(q), cellvertices(p)) == 1) {
                    PlaceHolder = cellvertices(pos - 1);
                    cellvertices(pos - 1) = cellvertices(p);
                    cellvertices(p) = PlaceHolder;

                    temp = values(pos - 1);
                    values(pos - 1) = values(p);
                    values(p) = temp;
                }
            }
        }
        if (ShareSide(cellvertices(pos - 2), cellvertices(pos - 3)) != 1) {
            if (ShareSide(cellvertices(pos - 3), cellvertices(pos - 1)) == 1 &&
                ShareSide(cellvertices(pos - N), cellvertices(pos - 2)) == 1) {
                PlaceHolder = cellvertices(pos - 2);
                cellvertices(pos - 2) = cellvertices(pos - 1);
                cellvertices(pos - 1) = PlaceHolder;

                temp = values(pos - 2);
                values(pos - 2) = values(pos - 1);
                values(pos - 1) = temp;
            }
        }
        if (ShareSide(cellvertices(pos - 1), cellvertices(pos - 2)) != 1) {
            if (ShareSide(cellvertices(pos - 3), cellvertices(pos - 1)) == 1 &&
                ShareSide(cellvertices(pos - 4), cellvertices(pos - 2)) == 1) {
                PlaceHolder = cellvertices(pos - 3);
                cellvertices(pos - 3) = cellvertices(pos - 2);
                cellvertices(pos - 2) = PlaceHolder;

                temp = values(pos - 3);
                values(pos - 3) = values(pos - 2);
                values(pos - 2) = temp;
            }
            if (ShareSide(cellvertices(pos - N + 1), cellvertices(pos - 3)) ==
                    1 &&
                ShareSide(cellvertices(pos - 1), cellvertices(pos - N + 1)) ==
                    1) {
                PlaceHolder = cellvertices(pos - 2);
                cellvertices(pos - 2) = cellvertices(pos - N + 1);
                cellvertices(pos - N + 1) = PlaceHolder;

                temp = values(pos - 2);
                values(pos - 2) = values(pos - N + 1);
                values(pos - N + 1) = temp;
            }
        }
        if (ShareSide(cellvertices(pos - N), cellvertices(pos - N + 1)) != 1) {
            if (ShareSide(cellvertices(pos - N), cellvertices(pos - 2)) == 1 &&
                ShareSide(cellvertices(pos - 1), cellvertices(pos - N + 1)) ==
                    1) {
                PlaceHolder = cellvertices(pos - 1);
                cellvertices(pos - 1) = cellvertices(pos - N);
                cellvertices(pos - N) = PlaceHolder;

                temp = values(pos - 1);
                values(pos - 1) = values(pos - N);
                values(pos - N) = temp;
            }
        }
    }

    // *    *    *   ESTABLISH TRIANGLE CONNECTIONS	   *    *    *

    for (p = pos - N + 2; p < pos; p++) {
        Tlist(0, nTris) = pos - N;
        Tlist(1, nTris) = p - 1;
        Tlist(2, nTris) = p;
        nTris++;
    }
}
//-------------------------------------------------------------------------------
/*inline void TRIM(DTMutableList<Point> &local_sol_pts, int &n_local_sol_pts, double isovalue,
				 IntArray &local_sol_tris, int &n_local_sol_tris,
				 DTMutableList<Point> &ns_pts, int &n_ns_pts, IntArray &ns_tris,
				 int &n_ns_tris, DTMutableList<Point> &ws_pts, int &n_ws_pts,
				 IntArray &ws_tris, int &n_ws_tris, DoubleArray &values,
				 DTMutableList<Point> &local_nws_pts, int &n_local_nws_pts,
				 DoubleArray &fluid_pad, DoubleArray &S, int &i, int &j, int &k,
				 int &newton_steps)
inline void TRIM(DTMutableList<DTPoint3D> &local_sol_pts, int &n_local_sol_pts, double isovalue,
				 DTMutableIntArray &local_sol_tris, int &n_local_sol_tris, 
				 DTMutableList<DTPoint3D> &ns_pts, int &n_ns_pts, DTMutableIntArray &ns_tris,
				 int &n_ns_tris, DTMutableList<DTPoint3D> &ws_pts, int &n_ws_pts, 
				 DTMutableIntArray &ws_tris, int &n_ws_tris, DTFloatArray &values,
				 DTMutableList<DTPoint3D> &local_nws_pts, int &n_local_nws_pts)
*/
//-------------------------------------------------------------------------------
inline void TRIM(DTMutableList<Point> &local_sol_pts, int &n_local_sol_pts,
                 double isovalue, IntArray &local_sol_tris,
                 int &n_local_sol_tris, DTMutableList<Point> &ns_pts,
                 int &n_ns_pts, IntArray &ns_tris, int &n_ns_tris,
                 DTMutableList<Point> &ws_pts, int &n_ws_pts, IntArray &ws_tris,
                 int &n_ws_tris, DoubleArray &values,
                 DTMutableList<Point> &local_nws_pts, int &n_local_nws_pts) {
    // Trim the local solid surface

    int map_ws;
    int map_ns;
    int p;
    int q;

    int a;
    int b;
    int c;
    Point A;
    Point B;
    Point C;
    Point D;
    Point E;
    Point P;

    bool all_to_ns = 1;
    // Check to see if all triangles in the cell are in ns_surface
    for (q = 0; q < n_local_sol_pts; q++) {
        if (values(q) < isovalue && all_to_ns == 1) {
            all_to_ns = 0;
        }
    }
    bool all_to_ws = 1;
    // Check to see if all triangles in the cell are in ws surface
    for (q = 0; q < n_local_sol_pts; q++) {
        if (values(q) > isovalue && all_to_ws == 1) {
            all_to_ws = 0;
        }
    }

    if (all_to_ws == 1) {
        map_ws = n_ws_pts;
        map_ws = 0;
        for (p = 0; p < n_local_sol_pts; p++) {
            ws_pts(n_ws_pts++) = local_sol_pts(p);
        }
        for (p = 0; p < n_local_sol_tris; p++) {
            ws_tris(0, n_ws_tris) = local_sol_tris(0, p) + map_ws;
            ws_tris(1, n_ws_tris) = local_sol_tris(1, p) + map_ws;
            ws_tris(2, n_ws_tris) = local_sol_tris(2, p) + map_ws;
            n_ws_tris++;
        }
    } else if (all_to_ns == 1) {
        map_ns = n_ns_pts;
        map_ns = 0;
        for (p = 0; p < n_local_sol_pts; p++) {
            ns_pts(n_ns_pts++) = local_sol_pts(p);
        }
        for (p = 0; p < n_local_sol_tris; p++) {
            ns_tris(0, n_ns_tris) = local_sol_tris(0, p) + map_ns;
            ns_tris(1, n_ns_tris) = local_sol_tris(1, p) + map_ns;
            ns_tris(2, n_ns_tris) = local_sol_tris(2, p) + map_ns;
            n_ns_tris++;
        }
    } else {
        // this section of surface contains a common line
        map_ns = n_ns_tris;
        map_ws = n_ws_tris;
        // Go through all triangles
        for (p = 0; p < n_local_sol_tris; p++) {
            a = local_sol_tris(0, p);
            b = local_sol_tris(1, p);
            c = local_sol_tris(2, p);
            A = local_sol_pts(a);
            B = local_sol_pts(b);
            C = local_sol_pts(c);
            if (values(a) > isovalue && values(b) > isovalue &&
                values(c) > isovalue) {
                // Triangle is in ns surface
                // Add points
                ns_pts(n_ns_pts++) = A;
                ns_pts(n_ns_pts++) = B;
                ns_pts(n_ns_pts++) = C;
                // Add triangles
                ns_tris(0, n_ns_tris) = n_ns_pts - 3;
                ns_tris(1, n_ns_tris) = n_ns_pts - 2;
                ns_tris(2, n_ns_tris) = n_ns_pts - 1;
                n_ns_tris++;
            } else if (values(a) < isovalue && values(b) < isovalue &&
                       values(c) < isovalue) {
                // Triangle is in ws surface
                // Add points
                ws_pts(n_ws_pts++) = A;
                ws_pts(n_ws_pts++) = B;
                ws_pts(n_ws_pts++) = C;
                // Add triangles
                ws_tris(0, n_ws_tris) = n_ws_pts - 3;
                ws_tris(1, n_ws_tris) = n_ws_pts - 2;
                ws_tris(2, n_ws_tris) = n_ws_pts - 1;
                n_ws_tris++;
            } else {
                // Triangle contains common line

                ////////////////////////////////////////
                ///////// FIND THE COMMON LINE /////////
                ////////////////////////////////////////

                if ((values(a) - isovalue) * (values(b) - isovalue) < 0) {
                    // compute common line vertex
                    //		P = A + (values(a) - isovalue)/(values(a)-values(b))*(B-A);
                    P.x = A.x + (values(a) - isovalue) /
                                    (values(a) - values(b)) * (B.x - A.x);
                    P.y = A.y + (values(a) - isovalue) /
                                    (values(a) - values(b)) * (B.y - A.y);
                    P.z = A.z + (values(a) - isovalue) /
                                    (values(a) - values(b)) * (B.z - A.z);

                    local_nws_pts(n_local_nws_pts++) = P;
/*					if ( n_local_nws_pts == 0 ){
						local_nws_pts(n_local_nws_pts++) = P;
					}
					else if ( P != local_nws_pts(n_local_nws_pts-1) ){
						local_nws_pts(n_local_nws_pts++) = P;
					}
*/				}
if ((values(b) - isovalue) * (values(c) - isovalue) < 0) {
    // compute common line vertex
    //					P = B + (values(b) - isovalue)/(values(b)-values(c))*(C-B);
    P.x = B.x + (values(b) - isovalue) / (values(b) - values(c)) * (C.x - B.x);
    P.y = B.y + (values(b) - isovalue) / (values(b) - values(c)) * (C.y - B.y);
    P.z = B.z + (values(b) - isovalue) / (values(b) - values(c)) * (C.z - B.z);

    local_nws_pts(n_local_nws_pts++) = P;
/*					if ( n_local_nws_pts == 0 ){
						local_nws_pts(n_local_nws_pts++) = P;
					}
					else if ( P != local_nws_pts(n_local_nws_pts-1) ){
						local_nws_pts(n_local_nws_pts++) = P;
					}
*/				}
if ((values(a) - isovalue) * (values(c) - isovalue) < 0) {
    // compute common line vertex
    //					P = A + (values(a) - isovalue)/(values(a)-values(c))*(C-A);
    P.x = A.x + (values(a) - isovalue) / (values(a) - values(c)) * (C.x - A.x);
    P.y = A.y + (values(a) - isovalue) / (values(a) - values(c)) * (C.y - A.y);
    P.z = A.z + (values(a) - isovalue) / (values(a) - values(c)) * (C.z - A.z);
    local_nws_pts(n_local_nws_pts++) = P;
}
// Store common line points as D and E
D = local_nws_pts(n_local_nws_pts - 2);
E = local_nws_pts(n_local_nws_pts - 1);
// Construct the new triangles
if ((values(a) - isovalue) * (values(b) - isovalue) < 0 &&
    (values(b) - isovalue) * (values(c) - isovalue) < 0) {
    if (values(b) > isovalue) {
        // Points
        ns_pts(n_ns_pts++) = B;
        ns_pts(n_ns_pts++) = D;
        ns_pts(n_ns_pts++) = E;
        // Triangles
        ns_tris(0, n_ns_tris) = n_ns_pts - 3; // B
        ns_tris(1, n_ns_tris) = n_ns_pts - 2; // D
        ns_tris(2, n_ns_tris) = n_ns_pts - 1; // E
        n_ns_tris++;
        // Points
        ws_pts(n_ws_pts++) = A;
        ws_pts(n_ws_pts++) = C;
        ws_pts(n_ws_pts++) = D;
        ws_pts(n_ws_pts++) = E;
        // Triangles (A,C,D),(C,D,E)
        ws_tris(0, n_ws_tris) = n_ws_pts - 4; // A
        ws_tris(1, n_ws_tris) = n_ws_pts - 3; // C
        ws_tris(2, n_ws_tris) = n_ws_pts - 2; // D
        n_ws_tris++;
        ws_tris(0, n_ws_tris) = n_ws_pts - 3; // C
        ws_tris(1, n_ws_tris) = n_ws_pts - 2; // D
        ws_tris(2, n_ws_tris) = n_ws_pts - 1; // E
        n_ws_tris++;
    } else {
        // Points
        ws_pts(n_ws_pts++) = B;
        ws_pts(n_ws_pts++) = D;
        ws_pts(n_ws_pts++) = E;
        // Triangles
        ws_tris(0, n_ws_tris) = n_ws_pts - 3; // B
        ws_tris(1, n_ws_tris) = n_ws_pts - 2; // D
        ws_tris(2, n_ws_tris) = n_ws_pts - 1; // E
        n_ws_tris++;
        // Points
        ns_pts(n_ns_pts++) = A;
        ns_pts(n_ns_pts++) = C;
        ns_pts(n_ns_pts++) = D;
        ns_pts(n_ns_pts++) = E;
        // Triangles (A,C,D),(C,D,E)
        ns_tris(0, n_ns_tris) = n_ns_pts - 4; // A
        ns_tris(1, n_ns_tris) = n_ns_pts - 3; // C
        ns_tris(2, n_ns_tris) = n_ns_pts - 2; // D
        n_ns_tris++;
        ns_tris(0, n_ns_tris) = n_ns_pts - 3; // C
        ns_tris(1, n_ns_tris) = n_ns_pts - 1; // E
        ns_tris(2, n_ns_tris) = n_ns_pts - 2; // D
        n_ns_tris++;
    }
} else if ((values(a) - isovalue) * (values(b) - isovalue) < 0 &&
           (values(a) - isovalue) * (values(c) - isovalue) < 0) {
    if (values(a) > isovalue) {
        // Points
        ns_pts(n_ns_pts++) = A;
        ns_pts(n_ns_pts++) = D;
        ns_pts(n_ns_pts++) = E;
        // Triangles
        ns_tris(0, n_ns_tris) = n_ns_pts - 3; // A
        ns_tris(1, n_ns_tris) = n_ns_pts - 2; // D
        ns_tris(2, n_ns_tris) = n_ns_pts - 1; // E
        n_ns_tris++;
        // Points
        ws_pts(n_ws_pts++) = B;
        ws_pts(n_ws_pts++) = C;
        ws_pts(n_ws_pts++) = D;
        ws_pts(n_ws_pts++) = E;
        // Triangles (B,C,D),(C,D,E)
        ws_tris(0, n_ws_tris) = n_ws_pts - 4; // B
        ws_tris(1, n_ws_tris) = n_ws_pts - 2; // D
        ws_tris(2, n_ws_tris) = n_ws_pts - 3; // C
        n_ws_tris++;
        ws_tris(0, n_ws_tris) = n_ws_pts - 3; // C
        ws_tris(1, n_ws_tris) = n_ws_pts - 2; // D
        ws_tris(2, n_ws_tris) = n_ws_pts - 1; // E
        n_ws_tris++;
    } else {
        // Points
        ws_pts(n_ws_pts++) = A;
        ws_pts(n_ws_pts++) = D;
        ws_pts(n_ws_pts++) = E;
        // Triangles
        ws_tris(0, n_ws_tris) = n_ws_pts - 3; // A
        ws_tris(1, n_ws_tris) = n_ws_pts - 2; // D
        ws_tris(2, n_ws_tris) = n_ws_pts - 1; // E
        n_ws_tris++;
        // Points
        ns_pts(n_ns_pts++) = B;
        ns_pts(n_ns_pts++) = C;
        ns_pts(n_ns_pts++) = D;
        ns_pts(n_ns_pts++) = E;
        // Triangles (B,C,D),(C,D,E)
        ns_tris(0, n_ns_tris) = n_ns_pts - 4; // B
        ns_tris(1, n_ns_tris) = n_ns_pts - 3; // C
        ns_tris(2, n_ns_tris) = n_ns_pts - 2; // D
        n_ns_tris++;
        ns_tris(0, n_ns_tris) = n_ns_pts - 3; // C
        ns_tris(1, n_ns_tris) = n_ns_pts - 1; // E
        ns_tris(2, n_ns_tris) = n_ns_pts - 2; // D
        n_ns_tris++;
    }
} else {
    if (values(c) > isovalue) {
        // Points
        ns_pts(n_ns_pts++) = C;
        ns_pts(n_ns_pts++) = D;
        ns_pts(n_ns_pts++) = E;
        // Triangles
        ns_tris(0, n_ns_tris) = n_ns_pts - 3; // C
        ns_tris(1, n_ns_tris) = n_ns_pts - 2; // D
        ns_tris(2, n_ns_tris) = n_ns_pts - 1; // E
        n_ns_tris++;
        // Points
        ws_pts(n_ws_pts++) = A;
        ws_pts(n_ws_pts++) = B;
        ws_pts(n_ws_pts++) = D;
        ws_pts(n_ws_pts++) = E;
        // Triangles (A,B,D),(A,D,E)
        ws_tris(0, n_ws_tris) = n_ws_pts - 4; // A
        ws_tris(1, n_ws_tris) = n_ws_pts - 3; // B
        ws_tris(2, n_ws_tris) = n_ws_pts - 2; // D
        n_ws_tris++;
        ws_tris(0, n_ws_tris) = n_ws_pts - 4; // A
        ws_tris(1, n_ws_tris) = n_ws_pts - 2; // D
        ws_tris(2, n_ws_tris) = n_ws_pts - 1; // E
        n_ws_tris++;
    } else {
        // Points
        ws_pts(n_ws_pts++) = C;
        ws_pts(n_ws_pts++) = D;
        ws_pts(n_ws_pts++) = E;
        // Triangles
        ws_tris(0, n_ws_tris) = n_ws_pts - 3; // C
        ws_tris(1, n_ws_tris) = n_ws_pts - 2; // D
        ws_tris(2, n_ws_tris) = n_ws_pts - 1; // E
        n_ws_tris++;
        // Points
        ns_pts(n_ns_pts++) = A;
        ns_pts(n_ns_pts++) = B;
        ns_pts(n_ns_pts++) = D;
        ns_pts(n_ns_pts++) = E;
        // Triangles (A,B,D),(A,D,E)
        ns_tris(0, n_ns_tris) = n_ns_pts - 4; // A
        ns_tris(1, n_ns_tris) = n_ns_pts - 3; // B
        ns_tris(2, n_ns_tris) = n_ns_pts - 2; // D
        n_ns_tris++;
        ns_tris(0, n_ns_tris) = n_ns_pts - 4; // A
        ns_tris(1, n_ns_tris) = n_ns_pts - 2; // D
        ns_tris(2, n_ns_tris) = n_ns_pts - 1; // E
        n_ns_tris++;
    }
}
            }
        }
    }
}
inline double geomavg_MarchingCubes(DoubleArray &A, double &v, int &i, int &j,
                                    int &k, DTMutableList<Point> &nw_pts,
                                    int &n_nw_pts, IntArray &nw_tris,
                                    int &n_nw_tris) {
    int N = 0; // n will be the number of vertices in this grid cell only
    Point P;
    Point pt;
    Point PlaceHolder;
    int m;
    int o;
    int p;

    // Go over each corner -- check to see if the corners are themselves vertices
    //1
    if (A(i, j, k) == v) {
        P.x = i;
        P.y = j;
        P.z = k;
        nw_pts(n_nw_pts++) = P;
        N++;
    }
    //2
    if (A(i + 1, j, k) == v) {
        P.x = i + 1;
        P.y = j;
        P.z = k;
        nw_pts(n_nw_pts++) = P;
        N++;
    }
    //3
    if (A(i + 1, j + 1, k) == v) {
        P.x = i + 1;
        P.y = j + 1;
        P.z = k;
        nw_pts(n_nw_pts++) = P;
        N++;
    }
    //4
    if (A(i, j + 1, k) == v) {
        P.x = i;
        P.y = j + 1;
        P.z = k;
        nw_pts(n_nw_pts++) = P;
        N++;
    }
    //5
    if (A(i, j, k + 1) == v) {
        P.x = i;
        P.y = j;
        P.z = k + 1;
        nw_pts(n_nw_pts++) = P;
        N++;
    }
    //6
    if (A(i + 1, j, k + 1) == v) {
        P.x = i + 1;
        P.y = j;
        P.z = k + 1;
        nw_pts(n_nw_pts++) = P;
        N++;
    }
    //7
    if (A(i + 1, j + 1, k + 1) == v) {
        P.x = i + 1;
        P.y = j + 1;
        P.z = k + 1;
        nw_pts(n_nw_pts++) = P;
        N++;
    }
    //8
    if (A(i, j + 1, k + 1) == v) {
        P.x = i;
        P.y = j + 1;
        P.z = k + 1;

        nw_pts(n_nw_pts++) = P;
        N++;
    }

    // Go through each side, compute P for sides of box spiraling up

    //	float val;
    if ((A(i, j, k) - v) * (A(i + 1, j, k) - v) < 0) {
        // If both points are in the fluid region
        if (A(i, j, k) != 0 && A(i + 1, j, k) != 0) {
            P.x = i + (A(i, j, k) - v) / (A(i, j, k) - A(i + 1, j, k));
            P.y = j;
            P.z = k;

            nw_pts(n_nw_pts++) = P;
            N++;
        }
    }

    if ((A(i + 1, j, k) - v) * (A(i + 1, j + 1, k) - v) < 0) {
        if (A(i + 1, j, k) != 0 && A(i + 1, j + 1, k) != 0) {
            P.x = i + 1;
            P.y = j +
                  (A(i + 1, j, k) - v) / (A(i + 1, j, k) - A(i + 1, j + 1, k));
            P.z = k;

            if (vertexcheck(P, N, n_nw_pts, nw_pts) ==
                1) { // P is a new vertex (not counted twice)
                nw_pts(n_nw_pts++) = P;
                N++;
            }
        }
    }

    if ((A(i + 1, j + 1, k) - v) * (A(i, j + 1, k) - v) < 0) {
        if (A(i + 1, j + 1, k) != 0 && A(i, j + 1, k) != 0) {
            P.x = i +
                  (A(i, j + 1, k) - v) / (A(i, j + 1, k) - A(i + 1, j + 1, k));
            P.y = j + 1;
            P.z = k;
            if (vertexcheck(P, N, n_nw_pts, nw_pts) ==
                1) { // P is a new vertex (not counted twice)
                nw_pts(n_nw_pts++) = P;
                N++;
            }
        }
    }

    //4
    if ((A(i, j + 1, k) - v) * (A(i, j, k) - v) < 0) {
        if (A(i, j + 1, k) != 0 && A(i, j, k) != 0) {
            P.x = i;
            P.y = j + (A(i, j, k) - v) / (A(i, j, k) - A(i, j + 1, k));
            P.z = k;
            if (vertexcheck(P, N, n_nw_pts, nw_pts) ==
                1) { // P is a new vertex (not counted twice)
                nw_pts(n_nw_pts++) = P;
                N++;
            }
        }
    }

    //5
    if ((A(i, j, k) - v) * (A(i, j, k + 1) - v) < 0) {
        if (A(i, j, k) != 0 && A(i, j, k + 1) != 0) {
            P.x = i;
            P.y = j;
            P.z = k + (A(i, j, k) - v) / (A(i, j, k) - A(i, j, k + 1));

            if (vertexcheck(P, N, n_nw_pts, nw_pts) ==
                1) { // P is a new vertex (not counted twice)
                nw_pts(n_nw_pts++) = P;
                N++;
            }
        }
    }

    //6
    if ((A(i + 1, j, k) - v) * (A(i + 1, j, k + 1) - v) < 0) {
        if (A(i + 1, j, k) != 0 && A(i + 1, j, k + 1) != 0) {
            P.x = i + 1;
            P.y = j;
            P.z = k +
                  (A(i + 1, j, k) - v) / (A(i + 1, j, k) - A(i + 1, j, k + 1));
            if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1) {
                nw_pts(n_nw_pts++) = P;
                N++;
            }
        }
    }

    //7
    if ((A(i + 1, j + 1, k) - v) * (A(i + 1, j + 1, k + 1) - v) < 0) {
        if (A(i + 1, j + 1, k) != 0 && A(i + 1, j + 1, k + 1) != 0) {
            P.x = i + 1;
            P.y = j + 1;
            P.z = k + (A(i + 1, j + 1, k) - v) /
                          (A(i + 1, j + 1, k) - A(i + 1, j + 1, k + 1));
            if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1) {
                nw_pts(n_nw_pts++) = P;
                N++;
            }
        }
    }

    //8
    if ((A(i, j + 1, k) - v) * (A(i, j + 1, k + 1) - v) < 0) {
        if (A(i, j + 1, k) != 0 && A(i, j + 1, k + 1) != 0) {
            P.x = i;
            P.y = j + 1;
            P.z = k +
                  (A(i, j + 1, k) - v) / (A(i, j + 1, k) - A(i, j + 1, k + 1));
            if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1) {
                nw_pts(n_nw_pts++) = P;
                N++;
            }
        }
    }

    //9
    if ((A(i, j, k + 1) - v) * (A(i + 1, j, k + 1) - v) < 0) {
        if (A(i, j, k + 1) != 0 && A(i + 1, j, k + 1) != 0) {
            P.x = i +
                  (A(i, j, k + 1) - v) / (A(i, j, k + 1) - A(i + 1, j, k + 1));
            P.y = j;
            P.z = k + 1;
            if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1) {
                nw_pts(n_nw_pts++) = P;
                N++;
            }
        }
    }

    //10
    if ((A(i + 1, j, k + 1) - v) * (A(i + 1, j + 1, k + 1) - v) < 0) {
        if (A(i + 1, j, k + 1) != 0 && A(i + 1, j + 1, k + 1) != 0) {
            P.x = i + 1;
            P.y = j + (A(i + 1, j, k + 1) - v) /
                          (A(i + 1, j, k + 1) - A(i + 1, j + 1, k + 1));
            P.z = k + 1;
            if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1) {
                nw_pts(n_nw_pts++) = P;
                N++;
            }
        }
    }

    //11
    if ((A(i + 1, j + 1, k + 1) - v) * (A(i, j + 1, k + 1) - v) < 0) {
        if (A(i + 1, j + 1, k + 1) != 0 && A(i, j + 1, k + 1) != 0) {
            P.x = i + (A(i, j + 1, k + 1) - v) /
                          (A(i, j + 1, k + 1) - A(i + 1, j + 1, k + 1));
            P.y = j + 1;
            P.z = k + 1;
            if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1) {
                nw_pts(n_nw_pts++) = P;
                N++;
            }
        }
    }

    //12
    if ((A(i, j + 1, k + 1) - v) * (A(i, j, k + 1) - v) < 0) {
        if (A(i, j + 1, k + 1) != 0 && A(i, j, k + 1) != 0) {
            P.x = i;
            P.y = j +
                  (A(i, j, k + 1) - v) / (A(i, j, k + 1) - A(i, j + 1, k + 1));
            P.z = k + 1;
            if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1) {
                nw_pts(n_nw_pts++) = P;
                N++;
            }
        }
    }

    // Assemble the triangles as long as points are found
    if (N > 0) {
        for (m = n_nw_pts - N; m < n_nw_pts - 2; m++) {
            for (o = m + 2; o < n_nw_pts - 1; o++) {
                if (ShareSide(nw_pts(m), nw_pts(o)) == 1) {
                    PlaceHolder = nw_pts(m + 1);
                    nw_pts(m + 1) = nw_pts(o);
                    nw_pts(o) = PlaceHolder;
                }
            }

            // make sure other neighbor of vertex 1 is in last spot
            if (m == n_nw_pts - N) {
                for (p = m + 2; p < n_nw_pts - 1; p++) {
                    if (ShareSide(nw_pts(m), nw_pts(p)) == 1) {
                        PlaceHolder = nw_pts(n_nw_pts - 1);
                        nw_pts(n_nw_pts - 1) = nw_pts(p);
                        nw_pts(p) = PlaceHolder;
                    }
                }
            }
            if (ShareSide(nw_pts(n_nw_pts - 2), nw_pts(n_nw_pts - 3)) != 1) {
                if (ShareSide(nw_pts(n_nw_pts - 3), nw_pts(n_nw_pts - 1)) ==
                        1 &&
                    ShareSide(nw_pts(n_nw_pts - N), nw_pts(n_nw_pts - 2)) ==
                        1) {
                    PlaceHolder = nw_pts(n_nw_pts - 2);
                    nw_pts(n_nw_pts - 2) = nw_pts(n_nw_pts - 1);
                    nw_pts(n_nw_pts - 1) = PlaceHolder;
                }
            }
            if (ShareSide(nw_pts(n_nw_pts - 1), nw_pts(n_nw_pts - 2)) != 1) {
                if (ShareSide(nw_pts(n_nw_pts - 3), nw_pts(n_nw_pts - 1)) ==
                        1 &&
                    ShareSide(nw_pts(n_nw_pts - 4), nw_pts(n_nw_pts - 2)) ==
                        1) {
                    PlaceHolder = nw_pts(n_nw_pts - 3);
                    nw_pts(n_nw_pts - 3) = nw_pts(n_nw_pts - 2);
                    nw_pts(n_nw_pts - 2) = PlaceHolder;
                }
                if (ShareSide(nw_pts(n_nw_pts - N + 1), nw_pts(n_nw_pts - 3)) ==
                        1 &&
                    ShareSide(nw_pts(n_nw_pts - 1), nw_pts(n_nw_pts - N + 1)) ==
                        1) {
                    PlaceHolder = nw_pts(n_nw_pts - 2);
                    nw_pts(n_nw_pts - 2) = nw_pts(n_nw_pts - N + 1);
                    nw_pts(n_nw_pts - N + 1) = PlaceHolder;
                }
            }
            if (ShareSide(nw_pts(n_nw_pts - N), nw_pts(n_nw_pts - N + 1)) !=
                1) {
                if (ShareSide(nw_pts(n_nw_pts - N), nw_pts(n_nw_pts - 2)) ==
                        1 &&
                    ShareSide(nw_pts(n_nw_pts - 1), nw_pts(n_nw_pts - N + 1)) ==
                        1) {
                    PlaceHolder = nw_pts(n_nw_pts - 1);
                    nw_pts(n_nw_pts - 1) = nw_pts(n_nw_pts - N);
                    nw_pts(n_nw_pts - N) = PlaceHolder;
                }
            }
        }

        // *    *    *   ESTABLISH TRIANGLE CONNECTIONS	   *    *    *

        for (p = n_nw_pts - N + 2; p < n_nw_pts; p++) {
            nw_tris(0, n_nw_tris) = n_nw_pts - N;
            nw_tris(1, n_nw_tris) = p - 1;
            nw_tris(2, n_nw_tris) = p;
            n_nw_tris++;
        }
    }
    // Compute the Interfacial Area
    double s1, s2, s3, s;
    Point pA, pB, pC;
    double area = 0.0;
    for (int r = 0; r < n_nw_tris; r++) {
        pA = nw_pts(nw_tris(0, r));
        pB = nw_pts(nw_tris(1, r));
        pC = nw_pts(nw_tris(2, r));
        // Compute length of sides (assume dx=dy=dz)
        s1 =
            sqrt((pA.x - pB.x) * (pA.x - pB.x) + (pA.y - pB.y) * (pA.y - pB.y) +
                 (pA.z - pB.z) * (pA.z - pB.z));
        s2 =
            sqrt((pA.x - pC.x) * (pA.x - pC.x) + (pA.y - pC.y) * (pA.y - pC.y) +
                 (pA.z - pC.z) * (pA.z - pC.z));
        s3 =
            sqrt((pB.x - pC.x) * (pB.x - pC.x) + (pB.y - pC.y) * (pB.y - pC.y) +
                 (pB.z - pC.z) * (pB.z - pC.z));
        s = 0.5 * (s1 + s2 + s3);
        area += sqrt(s * (s - s1) * (s - s2) * (s - s3));
    }
    return area;
}
//-------------------------------------------------------------------------------
inline void MC(DoubleArray &A, double &v, DoubleArray &solid, int &i, int &j,
               int &k, DTMutableList<Point> &nw_pts, int &n_nw_pts,
               IntArray &nw_tris, int &n_nw_tris) {
    int N = 0; // n will be the number of vertices in this grid cell only
    Point P;
    Point pt;
    DoubleArray TEST(3);

    Point PlaceHolder;
    int m;
    int o;
    int p;

    // Go over each corner -- check to see if the corners are themselves vertices
    //1
    if (A(i, j, k) == v) {
        P.x = i;
        P.y = j;
        P.z = k;
        nw_pts(n_nw_pts++) = P;
        N++;
    }
    //2
    if (A(i + 1, j, k) == v) {
        P.x = i + 1;
        P.y = j;
        P.z = k;
        nw_pts(n_nw_pts++) = P;
        N++;
    }
    //3
    if (A(i + 1, j + 1, k) == v) {
        P.x = i + 1;
        P.y = j + 1;
        P.z = k;
        nw_pts(n_nw_pts++) = P;
        N++;
    }
    //4
    if (A(i, j + 1, k) == v) {
        P.x = i;
        P.y = j + 1;
        P.z = k;
        nw_pts(n_nw_pts++) = P;
        N++;
    }
    //5
    if (A(i, j, k + 1) == v) {
        P.x = i;
        P.y = j;
        P.z = k + 1;
        nw_pts(n_nw_pts++) = P;
        N++;
    }
    //6
    if (A(i + 1, j, k + 1) == v) {
        P.x = i + 1;
        P.y = j;
        P.z = k + 1;
        nw_pts(n_nw_pts++) = P;
        N++;
    }
    //7
    if (A(i + 1, j + 1, k + 1) == v) {
        P.x = i + 1;
        P.y = j + 1;
        P.z = k + 1;
        nw_pts(n_nw_pts++) = P;
        N++;
    }
    //8
    if (A(i, j + 1, k + 1) == v) {
        P.x = i;
        P.y = j + 1;
        P.z = k + 1;

        nw_pts(n_nw_pts++) = P;
        N++;
    }

    // Go through each side, compute P for sides of box spiraling up

    //	float val;
    if ((A(i, j, k) - v) * (A(i + 1, j, k) - v) < 0) {
        // If both points are in the fluid region
        if (A(i, j, k) != 0 && A(i + 1, j, k) != 0) {
            P.x = i + (A(i, j, k) - v) / (A(i, j, k) - A(i + 1, j, k));
            P.y = j;
            P.z = k;
            // Evaluate the function S at the new point
            if (solid(i, j, k) * (1 - P.x + i) +
                    solid(i + 1, j, k) * (P.x - i) >
                0) {
                // This point is in the fluid region
                nw_pts(n_nw_pts++) = P;
                N++;
            }
        }
    }

    if ((A(i + 1, j, k) - v) * (A(i + 1, j + 1, k) - v) < 0) {
        if (A(i + 1, j, k) != 0 && A(i + 1, j + 1, k) != 0) {
            P.x = i + 1;
            P.y = j +
                  (A(i + 1, j, k) - v) / (A(i + 1, j, k) - A(i + 1, j + 1, k));
            P.z = k;
            // Evaluate the function S at the new point
            if (solid(i + 1, j, k) * (1 - P.y + j) +
                    solid(i + 1, j + 1, k) * (P.y - j) >
                0) {
                // This point is in the fluid region
                if (vertexcheck(P, N, n_nw_pts, nw_pts) ==
                    1) { // P is a new vertex (not counted twice)
                    nw_pts(n_nw_pts++) = P;
                    N++;
                }
            }
        }
    }

    if ((A(i + 1, j + 1, k) - v) * (A(i, j + 1, k) - v) < 0) {
        if (A(i + 1, j + 1, k) != 0 && A(i, j + 1, k) != 0) {
            P.x = i +
                  (A(i, j + 1, k) - v) / (A(i, j + 1, k) - A(i + 1, j + 1, k));
            P.y = j + 1;
            P.z = k;
            // Evaluate the function S at the new point
            if (solid(i, j + 1, k) * (1 - P.x + i) +
                    solid(i + 1, j + 1, k) * (P.x - i) >
                0) {
                // This point is in the fluid region
                if (vertexcheck(P, N, n_nw_pts, nw_pts) ==
                    1) { // P is a new vertex (not counted twice)
                    nw_pts(n_nw_pts++) = P;
                    N++;
                }
            }
        }
    }

    //4
    if ((A(i, j + 1, k) - v) * (A(i, j, k) - v) < 0) {
        if (A(i, j + 1, k) != 0 && A(i, j, k) != 0) {
            P.x = i;
            P.y = j + (A(i, j, k) - v) / (A(i, j, k) - A(i, j + 1, k));
            P.z = k;
            // Evaluate the function S at the new point
            if (solid(i, j, k) * (1 - P.y + j) +
                    solid(i, j + 1, k) * (P.y - j) >
                0) {
                // This point is in the fluid region
                if (vertexcheck(P, N, n_nw_pts, nw_pts) ==
                    1) { // P is a new vertex (not counted twice)
                    nw_pts(n_nw_pts++) = P;
                    N++;
                }
            }
        }
    }

    //5
    if ((A(i, j, k) - v) * (A(i, j, k + 1) - v) < 0) {
        if (A(i, j, k) != 0 && A(i, j, k + 1) != 0) {
            P.x = i;
            P.y = j;
            P.z = k + (A(i, j, k) - v) / (A(i, j, k) - A(i, j, k + 1));
            // Evaluate the function S at the new point
            if (solid(i, j, k) * (1 - P.z + k) +
                    solid(i, j, k + 1) * (P.z - k) >
                0) {
                // This point is in the fluid region
                if (vertexcheck(P, N, n_nw_pts, nw_pts) ==
                    1) { // P is a new vertex (not counted twice)
                    nw_pts(n_nw_pts++) = P;
                    N++;
                }
            }
        }
    }

    //6
    if ((A(i + 1, j, k) - v) * (A(i + 1, j, k + 1) - v) < 0) {
        if (A(i + 1, j, k) != 0 && A(i + 1, j, k + 1) != 0) {
            P.x = i + 1;
            P.y = j;
            P.z = k +
                  (A(i + 1, j, k) - v) / (A(i + 1, j, k) - A(i + 1, j, k + 1));
            // Evaluate the function S at the new point
            if (solid(i + 1, j, k) * (1 - P.z + k) +
                    solid(i + 1, j, k + 1) * (P.z - k) >
                0) {
                // This point is in the fluid region
                if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1) {
                    nw_pts(n_nw_pts++) = P;
                    N++;
                }
            }
        }
    }

    //7
    if ((A(i + 1, j + 1, k) - v) * (A(i + 1, j + 1, k + 1) - v) < 0) {
        if (A(i + 1, j + 1, k) != 0 && A(i + 1, j + 1, k + 1) != 0) {
            P.x = i + 1;
            P.y = j + 1;
            P.z = k + (A(i + 1, j + 1, k) - v) /
                          (A(i + 1, j + 1, k) - A(i + 1, j + 1, k + 1));
            // Evaluate the function S at the new point
            if (solid(i + 1, j + 1, k) * (1 - P.z + k) +
                    solid(i + 1, j + 1, k + 1) * (P.z - k) >
                0) {
                // This point is in the fluid region
                if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1) {
                    nw_pts(n_nw_pts++) = P;
                    N++;
                }
            }
        }
    }

    //8
    if ((A(i, j + 1, k) - v) * (A(i, j + 1, k + 1) - v) < 0) {
        if (A(i, j + 1, k) != 0 && A(i, j + 1, k + 1) != 0) {
            P.x = i;
            P.y = j + 1;
            P.z = k +
                  (A(i, j + 1, k) - v) / (A(i, j + 1, k) - A(i, j + 1, k + 1));
            // Evaluate the function S at the new point
            if (solid(i, j + 1, k) * (1 - P.z + k) +
                    solid(i, j + 1, k + 1) * (P.z - k) >
                0) {
                // This point is in the fluid region
                if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1) {
                    nw_pts(n_nw_pts++) = P;
                    N++;
                }
            }
        }
    }

    //9
    if ((A(i, j, k + 1) - v) * (A(i + 1, j, k + 1) - v) < 0) {
        if (A(i, j, k + 1) != 0 && A(i + 1, j, k + 1) != 0) {
            P.x = i +
                  (A(i, j, k + 1) - v) / (A(i, j, k + 1) - A(i + 1, j, k + 1));
            P.y = j;
            P.z = k + 1;
            // Evaluate the function S at the new point
            if (solid(i, j, k + 1) * (1 - P.x + i) +
                    solid(i + 1, j, k + 1) * (P.x - i) >
                0) {
                // This point is in the fluid region
                if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1) {
                    nw_pts(n_nw_pts++) = P;
                    N++;
                }
            }
        }
    }

    //10
    if ((A(i + 1, j, k + 1) - v) * (A(i + 1, j + 1, k + 1) - v) < 0) {
        if (A(i + 1, j, k + 1) != 0 && A(i + 1, j + 1, k + 1) != 0) {
            P.x = i + 1;
            P.y = j + (A(i + 1, j, k + 1) - v) /
                          (A(i + 1, j, k + 1) - A(i + 1, j + 1, k + 1));
            P.z = k + 1;
            // Evaluate the function S at the new point
            if (solid(i + 1, j, k + 1) * (1 - P.y + j) +
                    solid(i + 1, j + 1, k + 1) * (P.y - j) >
                0) {
                // This point is in the fluid region
                if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1) {
                    nw_pts(n_nw_pts++) = P;
                    N++;
                }
            }
        }
    }

    //11
    if ((A(i + 1, j + 1, k + 1) - v) * (A(i, j + 1, k + 1) - v) < 0) {
        if (A(i + 1, j + 1, k + 1) != 0 && A(i, j + 1, k + 1) != 0) {
            P.x = i + (A(i, j + 1, k + 1) - v) /
                          (A(i, j + 1, k + 1) - A(i + 1, j + 1, k + 1));
            P.y = j + 1;
            P.z = k + 1;
            // Evaluate the function S at the new point
            if (solid(i, j + 1, k + 1) * (1 - P.x + i) +
                    solid(i + 1, j + 1, k + 1) * (P.x - i) >
                0) {
                // This point is in the fluid region
                if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1) {
                    nw_pts(n_nw_pts++) = P;
                    N++;
                }
            }
        }
    }

    //12
    if ((A(i, j + 1, k + 1) - v) * (A(i, j, k + 1) - v) < 0) {
        if (A(i, j + 1, k + 1) != 0 && A(i, j, k + 1) != 0) {
            P.x = i;
            P.y = j +
                  (A(i, j, k + 1) - v) / (A(i, j, k + 1) - A(i, j + 1, k + 1));
            P.z = k + 1;
            // Evaluate the function S at the new point
            if (solid(i, j, k + 1) * (1 - P.y + j) +
                    solid(i, j + 1, k + 1) * (P.y - j) >
                0) {
                // This point is in the fluid region
                if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1) {
                    nw_pts(n_nw_pts++) = P;
                    N++;
                }
            }
        }
    }

    // sort vertices so that they are connected to "neighbors"

    // DTMatlabDataFile("/tmp/Dump.mat",DTFile::NewReadWrite);

    // n_nw_pts = number of vertices in total (location n_nw_pts is first unfilled position)
    // n = number of vertices in this grid cell

    // Assemble the triangles as long as points are found
    if (N > 0) {
        for (m = n_nw_pts - N; m < n_nw_pts - 2; m++) {
            for (o = m + 2; o < n_nw_pts - 1; o++) {
                if (ShareSide(nw_pts(m), nw_pts(o)) == 1) {
                    PlaceHolder = nw_pts(m + 1);
                    nw_pts(m + 1) = nw_pts(o);
                    nw_pts(o) = PlaceHolder;
                }
            }

            // make sure other neighbor of vertex 1 is in last spot
            if (m == n_nw_pts - N) {
                for (p = m + 2; p < n_nw_pts - 1; p++) {
                    if (ShareSide(nw_pts(m), nw_pts(p)) == 1) {
                        PlaceHolder = nw_pts(n_nw_pts - 1);
                        nw_pts(n_nw_pts - 1) = nw_pts(p);
                        nw_pts(p) = PlaceHolder;
                    }
                }
            }
            if (ShareSide(nw_pts(n_nw_pts - 2), nw_pts(n_nw_pts - 3)) != 1) {
                if (ShareSide(nw_pts(n_nw_pts - 3), nw_pts(n_nw_pts - 1)) ==
                        1 &&
                    ShareSide(nw_pts(n_nw_pts - N), nw_pts(n_nw_pts - 2)) ==
                        1) {
                    PlaceHolder = nw_pts(n_nw_pts - 2);
                    nw_pts(n_nw_pts - 2) = nw_pts(n_nw_pts - 1);
                    nw_pts(n_nw_pts - 1) = PlaceHolder;
                }
            }
            if (ShareSide(nw_pts(n_nw_pts - 1), nw_pts(n_nw_pts - 2)) != 1) {
                if (ShareSide(nw_pts(n_nw_pts - 3), nw_pts(n_nw_pts - 1)) ==
                        1 &&
                    ShareSide(nw_pts(n_nw_pts - 4), nw_pts(n_nw_pts - 2)) ==
                        1) {
                    PlaceHolder = nw_pts(n_nw_pts - 3);
                    nw_pts(n_nw_pts - 3) = nw_pts(n_nw_pts - 2);
                    nw_pts(n_nw_pts - 2) = PlaceHolder;
                }
                if (ShareSide(nw_pts(n_nw_pts - N + 1), nw_pts(n_nw_pts - 3)) ==
                        1 &&
                    ShareSide(nw_pts(n_nw_pts - 1), nw_pts(n_nw_pts - N + 1)) ==
                        1) {
                    PlaceHolder = nw_pts(n_nw_pts - 2);
                    nw_pts(n_nw_pts - 2) = nw_pts(n_nw_pts - N + 1);
                    nw_pts(n_nw_pts - N + 1) = PlaceHolder;
                }
            }
            if (ShareSide(nw_pts(n_nw_pts - N), nw_pts(n_nw_pts - N + 1)) !=
                1) {
                if (ShareSide(nw_pts(n_nw_pts - N), nw_pts(n_nw_pts - 2)) ==
                        1 &&
                    ShareSide(nw_pts(n_nw_pts - 1), nw_pts(n_nw_pts - N + 1)) ==
                        1) {
                    PlaceHolder = nw_pts(n_nw_pts - 1);
                    nw_pts(n_nw_pts - 1) = nw_pts(n_nw_pts - N);
                    nw_pts(n_nw_pts - N) = PlaceHolder;
                }
            }
        }

        // *    *    *   ESTABLISH TRIANGLE CONNECTIONS	   *    *    *

        for (p = n_nw_pts - N + 2; p < n_nw_pts; p++) {
            nw_tris(0, n_nw_tris) = n_nw_pts - N;
            nw_tris(1, n_nw_tris) = p - 1;
            nw_tris(2, n_nw_tris) = p;
            n_nw_tris++;
        }
    }
}
//-------------------------------------------------------------------------------
inline void EDGE(DoubleArray &A, double &v, DoubleArray &solid, int &i, int &j,
                 int &k, int &m, int &n, int &o, DTMutableList<Point> &nw_pts,
                 int &n_nw_pts, IntArray &nw_tris, int &n_nw_tris,
                 DTMutableList<Point> &local_nws_pts, int &n_local_nws_pts) {

    NULL_USE(m);
    NULL_USE(n);
    NULL_USE(o);

    // FIND THE POINTS ON THE nw SURFACE THAT ARE ON THE EDGE (COMMON LINE WITH SOLID PHASE)
    // function A is the fluid data padded (so that it has values inside the solid phase)

    int N = 0; // n will be the number of vertices in this grid cell only
    Point P;

    Point temp;

    int p;
    int q;
    int r;

    // Add common line points to nw_pts
    for (p = 0; p < n_local_nws_pts; p++) {
        nw_pts(n_nw_pts++) = local_nws_pts(p);
    }

    // Go over each corner -- check to see if the corners are themselves vertices
    //1
    if (A(i, j, k) == v) {
        P.x = i;
        P.y = j;
        P.z = k;
        nw_pts(n_nw_pts++) = P;
        N++;
    }
    //2
    if (A(i + 1, j, k) == v) {
        P.x = i + 1;
        P.y = j;
        P.z = k;
        nw_pts(n_nw_pts++) = P;
        N++;
    }
    //3
    if (A(i + 1, j + 1, k) == v) {
        P.x = i + 1;
        P.y = j + 1;
        P.z = k;
        nw_pts(n_nw_pts++) = P;
        N++;
    }
    //4
    if (A(i, j + 1, k) == v) {
        P.x = i;
        P.y = j + 1;
        P.z = k;
        nw_pts(n_nw_pts++) = P;
        N++;
    }
    //5
    if (A(i, j, k + 1) == v) {
        P.x = i;
        P.y = j;
        P.z = k + 1;
        nw_pts(n_nw_pts++) = P;
        N++;
    }
    //6
    if (A(i + 1, j, k + 1) == v) {
        P.x = i + 1;
        P.y = j;
        P.z = k + 1;
        nw_pts(n_nw_pts++) = P;
        N++;
    }
    //7
    if (A(i + 1, j + 1, k + 1) == v) {
        P.x = i + 1;
        P.y = j + 1;
        P.z = k + 1;
        nw_pts(n_nw_pts++) = P;
        N++;
    }
    //8
    if (A(i, j + 1, k + 1) == v) {
        P.x = i;
        P.y = j + 1;
        P.z = k + 1;

        nw_pts(n_nw_pts++) = P;
        N++;
    }

    // Go through each side, compute P for sides of box spiraling up
    //	float val;
    Point pt;

    // 1
    if ((A(i, j, k) - v) * (A(i + 1, j, k) - v) < 0) {
        // If both points are in the fluid region
        if (A(i, j, k) != 0 && A(i + 1, j, k) != 0) {
            P.x = i + (A(i, j, k) - v) / (A(i, j, k) - A(i + 1, j, k));
            P.y = j;
            P.z = k;
            // Evaluate the function S at the new point
            if (solid(i, j, k) * (1 - P.x + i) +
                    solid(i + 1, j, k) * (P.x - i) >
                0) {
                // This point is in the fluid region
                nw_pts(n_nw_pts++) = P;
                N++;
            }
        }
    }

    // 2
    if ((A(i + 1, j, k) - v) * (A(i + 1, j + 1, k) - v) < 0) {
        if (A(i + 1, j, k) != 0 && A(i + 1, j + 1, k) != 0) {
            P.x = i + 1;
            P.y = j +
                  (A(i + 1, j, k) - v) / (A(i + 1, j, k) - A(i + 1, j + 1, k));
            P.z = k;
            // Evaluate the function S at the new point
            if (solid(i + 1, j, k) * (1 - P.y + j) +
                    solid(i + 1, j + 1, k) * (P.y - j) >
                0) {
                // This point is in the fluid region
                if (vertexcheck(P, N, n_nw_pts, nw_pts) ==
                    1) { // P is a new vertex (not counted twice)
                    nw_pts(n_nw_pts++) = P;
                    N++;
                }
            }
        }
    }

    //3
    if ((A(i + 1, j + 1, k) - v) * (A(i, j + 1, k) - v) < 0) {
        if (A(i + 1, j + 1, k) != 0 && A(i, j + 1, k) != 0) {
            P.x = i +
                  (A(i, j + 1, k) - v) / (A(i, j + 1, k) - A(i + 1, j + 1, k));
            P.y = j + 1;
            P.z = k;
            // Evaluate the function S at the new point
            if (solid(i, j + 1, k) * (1 - P.x + i) +
                    solid(i + 1, j + 1, k) * (P.x - i) >
                0) {
                // This point is in the fluid region
                if (vertexcheck(P, N, n_nw_pts, nw_pts) ==
                    1) { // P is a new vertex (not counted twice)
                    nw_pts(n_nw_pts++) = P;
                    N++;
                }
            }
        }
    }

    //4
    if ((A(i, j + 1, k) - v) * (A(i, j, k) - v) < 0) {
        if (A(i, j + 1, k) != 0 && A(i, j, k) != 0) {
            P.x = i;
            P.y = j + (A(i, j, k) - v) / (A(i, j, k) - A(i, j + 1, k));
            P.z = k;
            // Evaluate the function S at the new point
            if (solid(i, j, k) * (1 - P.y + j) +
                    solid(i, j + 1, k) * (P.y - j) >
                0) {
                // This point is in the fluid region
                if (vertexcheck(P, N, n_nw_pts, nw_pts) ==
                    1) { // P is a new vertex (not counted twice)
                    nw_pts(n_nw_pts++) = P;
                    N++;
                }
            }
        }
    }
    //5
    if ((A(i, j, k) - v) * (A(i, j, k + 1) - v) < 0) {
        if (A(i, j, k) != 0 && A(i, j, k + 1) != 0) {
            P.x = i;
            P.y = j;
            P.z = k + (A(i, j, k) - v) / (A(i, j, k) - A(i, j, k + 1));
            // Evaluate the function S at the new point
            if (solid(i, j, k) * (1 - P.z + k) +
                    solid(i, j, k + 1) * (P.z - k) >
                0) {
                // This point is in the fluid region
                if (vertexcheck(P, N, n_nw_pts, nw_pts) ==
                    1) { // P is a new vertex (not counted twice)
                    nw_pts(n_nw_pts++) = P;
                    N++;
                }
            }
        }
    }

    //6
    if ((A(i + 1, j, k) - v) * (A(i + 1, j, k + 1) - v) < 0) {
        if (A(i + 1, j, k) != 0 && A(i + 1, j, k + 1) != 0) {
            P.x = i + 1;
            P.y = j;
            P.z = k +
                  (A(i + 1, j, k) - v) / (A(i + 1, j, k) - A(i + 1, j, k + 1));
            // Evaluate the function S at the new point
            if (solid(i + 1, j, k) * (1 - P.z + k) +
                    solid(i + 1, j, k + 1) * (P.z - k) >
                0) {
                // This point is in the fluid region
                if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1) {
                    nw_pts(n_nw_pts++) = P;
                    N++;
                }
            }
        }
    }
    //7
    if ((A(i + 1, j + 1, k) - v) * (A(i + 1, j + 1, k + 1) - v) < 0) {
        if (A(i + 1, j + 1, k) != 0 && A(i + 1, j + 1, k + 1) != 0) {
            P.x = i + 1;
            P.y = j + 1;
            P.z = k + (A(i + 1, j + 1, k) - v) /
                          (A(i + 1, j + 1, k) - A(i + 1, j + 1, k + 1));
            // Evaluate the function S at the new point
            if (solid(i + 1, j + 1, k) * (1 - P.z + k) +
                    solid(i + 1, j + 1, k + 1) * (P.z - k) >
                0) {
                // This point is in the fluid region
                if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1) {
                    nw_pts(n_nw_pts++) = P;
                    N++;
                }
            }
        }
    }

    //8
    if ((A(i, j + 1, k) - v) * (A(i, j + 1, k + 1) - v) < 0) {
        if (A(i, j + 1, k) != 0 && A(i, j + 1, k + 1) != 0) {
            P.x = i;
            P.y = j + 1;
            P.z = k +
                  (A(i, j + 1, k) - v) / (A(i, j + 1, k) - A(i, j + 1, k + 1));
            // Evaluate the function S at the new point
            if (solid(i, j + 1, k) * (1 - P.z + k) +
                    solid(i, j + 1, k + 1) * (P.z - k) >
                0) {
                // This point is in the fluid region
                if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1) {
                    nw_pts(n_nw_pts++) = P;
                    N++;
                }
            }
        }
    }

    //9
    if ((A(i, j, k + 1) - v) * (A(i + 1, j, k + 1) - v) < 0) {
        if (A(i, j, k + 1) != 0 && A(i + 1, j, k + 1) != 0) {
            P.x = i +
                  (A(i, j, k + 1) - v) / (A(i, j, k + 1) - A(i + 1, j, k + 1));
            P.y = j;
            P.z = k + 1;
            // Evaluate the function S at the new point
            if (solid(i, j, k + 1) * (1 - P.x + i) +
                    solid(i + 1, j, k + 1) * (P.x - i) >
                0) {
                // This point is in the fluid region
                if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1) {
                    nw_pts(n_nw_pts++) = P;
                    N++;
                }
            }
        }
    }

    //10
    if ((A(i + 1, j, k + 1) - v) * (A(i + 1, j + 1, k + 1) - v) < 0) {
        if (A(i + 1, j, k + 1) != 0 && A(i + 1, j + 1, k + 1) != 0) {
            P.x = i + 1;
            P.y = j + (A(i + 1, j, k + 1) - v) /
                          (A(i + 1, j, k + 1) - A(i + 1, j + 1, k + 1));
            P.z = k + 1;
            // Evaluate the function S at the new point
            if (solid(i + 1, j, k + 1) * (1 - P.y + j) +
                    solid(i + 1, j + 1, k + 1) * (P.y - j) >
                0) {
                // This point is in the fluid region
                if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1) {
                    nw_pts(n_nw_pts++) = P;
                    N++;
                }
            }
        }
    }

    //11
    if ((A(i + 1, j + 1, k + 1) - v) * (A(i, j + 1, k + 1) - v) < 0) {
        if (A(i + 1, j + 1, k + 1) != 0 && A(i, j + 1, k + 1) != 0) {
            P.x = i + (A(i, j + 1, k + 1) - v) /
                          (A(i, j + 1, k + 1) - A(i + 1, j + 1, k + 1));
            P.y = j + 1;
            P.z = k + 1;
            // Evaluate the function S at the new point
            if (solid(i, j + 1, k + 1) * (1 - P.x + i) +
                    solid(i + 1, j + 1, k + 1) * (P.x - i) >
                0) {
                // This point is in the fluid region
                if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1) {
                    nw_pts(n_nw_pts++) = P;
                    N++;
                }
            }
        }
    }

    //12
    if ((A(i, j + 1, k + 1) - v) * (A(i, j, k + 1) - v) < 0) {
        if (A(i, j + 1, k + 1) != 0 && A(i, j, k + 1) != 0) {
            P.x = i;
            P.y = j +
                  (A(i, j, k + 1) - v) / (A(i, j, k + 1) - A(i, j + 1, k + 1));
            P.z = k + 1;
            // Evaluate the function S at the new point
            if (solid(i, j, k + 1) * (1 - P.y + j) +
                    solid(i, j + 1, k + 1) * (P.y - j) >
                0) {
                // This point is in the fluid region
                if (vertexcheck(P, N, n_nw_pts, nw_pts) == 1) {
                    nw_pts(n_nw_pts++) = P;
                    N++;
                }
            }
        }
    }

    ////////////////////////////////////////////////
    ////// SORT VERTICES SO THAT THEY CONNECT //////
    //////        TO ALL "NEIGHBORS"          //////
    ////////////////////////////////////////////////

    //  First common line point should connect to last MC point
    for (q = n_nw_pts - N; q < n_nw_pts - 1; q++) {
        if (ShareSide(nw_pts(n_nw_pts - N - n_local_nws_pts),
                      nw_pts(n_nw_pts - 1)) == 0 &&
            ShareSide(nw_pts(n_nw_pts - N - n_local_nws_pts), nw_pts(q)) == 1) {
            // Switch point q with point n_nw_pts-N-n_local_nws_pts
            temp = nw_pts(n_nw_pts - 1);
            nw_pts(n_nw_pts - 1) = nw_pts(q);
            nw_pts(q) = temp;
        }
    }
    //  Last common line point should connect to first MC point
    for (q = n_nw_pts - N + 1; q < n_nw_pts - 1; q++) {
        if (ShareSide(nw_pts(n_nw_pts - N - 1), nw_pts(n_nw_pts - N)) == 0 &&
            ShareSide(nw_pts(n_nw_pts - N - 1), nw_pts(q)) == 1) {
            // Switch these points
            temp = nw_pts(n_nw_pts - N);
            nw_pts(n_nw_pts - N) = nw_pts(q);
            nw_pts(q) = temp;
        }
    }
    // All MC points should connect to their neighbors
    for (q = n_nw_pts - N; q < n_nw_pts - 2; q++) {
        if (ShareSide(nw_pts(q), nw_pts(q + 1)) == 0) {
            for (r = q + 2; r < n_nw_pts - 1; r++) {
                if (ShareSide(nw_pts(q), nw_pts(q + 1)) == 0 &&
                    ShareSide(nw_pts(q), nw_pts(r)) == 1) {
                    // Switch r and q+1
                    temp = nw_pts(q + 1);
                    nw_pts(q + 1) = nw_pts(r);
                    nw_pts(r) = temp;
                }
            }
        }
    }

    // *    *    *   ESTABLISH TRIANGLE CONNECTIONS	   *    *    *
    for (p = n_nw_pts - N - n_local_nws_pts; p < n_nw_pts - 2; p++) {
        nw_tris(0, n_nw_tris) = n_nw_pts - 1;
        nw_tris(1, n_nw_tris) = p;
        nw_tris(2, n_nw_tris) = p + 1;
        n_nw_tris++;
    }

    //	for (p=n_nw_pts-N-n_local_nws_pts+2; p<n_nw_pts; p++){
    //		nw_tris(0,n_nw_tris) = n_nw_pts-N-n_local_nws_pts;
    //		nw_tris(1,n_nw_tris) = p-1;
    //		nw_tris(2,n_nw_tris) = p;
    //        n_nw_tris++;
    //	}
}
//--------------------------------------------------------------------------------------------------------

inline bool PointsEqual(const Point &A, const Point &B) {
    bool equal = false;

    //	double value;
    //	value = (A.x-B.x)*(A.x-B.x) + (A.y-B.y)*(A.y-B.y) + (A.z-B.z)*(A.z-B.z);
    //	if (value < 1e-14) equal = true;

    if (A.x == B.x && A.y == B.y && A.z == B.z)
        equal = true;

    return equal;
}
//--------------------------------------------------------------------------------------------------------
inline void ComputeAreasPMMC(IntArray &cubeList, int start, int finish,
                             DoubleArray &F, DoubleArray &S, double vF,
                             double vS, double &blob_volume, double &ans,
                             double &aws, double &awn, double &lwns, int Nx,
                             int Ny, int Nz) {
    /* ****************************************************************
	 VARIABLES FOR THE PMMC ALGORITHM
	 ****************************************************************** */
    awn = aws = ans = lwns = 0.0;

    //	bool add=1;			// Set to false if any corners contain nw-phase ( F > vF)
    int cube[8][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 1},
                      {1, 0, 1}, {0, 1, 1}, {1, 1, 1}}; // cube corners
    //	int count_in=0,count_out=0;
    //	int nodx,nody,nodz;
    // initialize lists for vertices for surfaces, common line
    DTMutableList<Point> nw_pts(20);
    DTMutableList<Point> ns_pts(20);
    DTMutableList<Point> ws_pts(20);
    DTMutableList<Point> nws_pts(20);
    // initialize triangle lists for surfaces
    IntArray nw_tris(3, 20);
    IntArray ns_tris(3, 20);
    IntArray ws_tris(3, 20);
    // initialize list for line segments
    IntArray nws_seg(2, 20);

    DTMutableList<Point> tmp(20);
    //	IntArray store;

    int i, j, k, p, q, r;
    int n_nw_pts = 0, n_ns_pts = 0, n_ws_pts = 0, n_nws_pts = 0, map = 0;
    int n_nw_tris = 0, n_ns_tris = 0, n_ws_tris = 0, n_nws_seg = 0;

    double s, s1, s2, s3; // Triangle sides (lengths)
    Point A, B, C, P;
    //	double area;

    // Initialize arrays for local solid surface
    DTMutableList<Point> local_sol_pts(20);
    int n_local_sol_pts = 0;
    IntArray local_sol_tris(3, 18);
    int n_local_sol_tris;
    DoubleArray values(20);
    DTMutableList<Point> local_nws_pts(20);
    int n_local_nws_pts;

    int n_nw_tris_beg, n_ns_tris_beg, n_ws_tris_beg;
    int c;
    //	int newton_steps = 0;
    //	double blob_volume;

    /* ****************************************************************
	 RUN PMMC ON EACH BLOB
	 ****************************************************************** */
    //	printf("Running the PMMC Algorithm \n");
    //	printf("The number of blobs is %i \n",nblobs);

    // Store beginning points for surfaces for blob p
    n_nw_tris_beg = n_nw_tris;
    n_ns_tris_beg = n_ns_tris;
    n_ws_tris_beg = n_ws_tris;
    //	n_nws_seg_beg = n_nws_seg;
    // Loop over all cubes
    blob_volume = 0; // Initialize the volume for blob a to zero
    for (c = start; c < finish; c++) {
        // Get cube from the list
        i = cubeList(0, c);
        j = cubeList(1, c);
        k = cubeList(2, c);

        //	printf("somewhere inside %i \n",c);

        /*		// Compute the volume
		for (p=0;p<8;p++){
			if ( indicator(i+cube[p][0],j+cube[p][1],k+cube[p][2]) > -1){
				blob_volume += 0.125;
			}
		}
*/
        for (p = 0; p < 8; p++) {
            if (F(i + cube[p][0], j + cube[p][1], k + cube[p][2]) > 0 &&
                S(i + cube[p][0], j + cube[p][1], k + cube[p][2]) > 0) {
                blob_volume += 0.125;
            }
        }

        // Run PMMC
        n_local_sol_tris = 0;
        n_local_sol_pts = 0;
        n_local_nws_pts = 0;

        // if there is a solid phase interface in the grid cell
        if (Interface(S, vS, i, j, k) == 1) {

            /////////////////////////////////////////
            /// CONSTRUCT THE LOCAL SOLID SURFACE ///
            /////////////////////////////////////////

            // find the local solid surface
            SOL_SURF(S, 0, F, vF, i, j, k, Nx, Ny, Nz, local_sol_pts,
                     n_local_sol_pts, local_sol_tris, n_local_sol_tris, values);

            /////////////////////////////////////////
            //////// TRIM THE SOLID SURFACE /////////
            /////////////////////////////////////////
            /*			TRIM(local_sol_pts, n_local_sol_pts, fluid_isovalue,local_sol_tris, n_local_sol_tris,
									 ns_pts, n_ns_pts, ns_tris, n_ns_tris, ws_pts, n_ws_pts,
									 ws_tris, n_ws_tris, values, local_nws_pts, n_local_nws_pts,
									 Phase, SignDist, i, j, k, newton_steps);
*/
            TRIM(local_sol_pts, n_local_sol_pts, vF, local_sol_tris,
                 n_local_sol_tris, ns_pts, n_ns_pts, ns_tris, n_ns_tris, ws_pts,
                 n_ws_pts, ws_tris, n_ws_tris, values, local_nws_pts,
                 n_local_nws_pts);

            /////////////////////////////////////////
            //////// WRITE COMMON LINE POINTS ///////
            ////////      TO MAIN ARRAYS      ///////
            /////////////////////////////////////////
            map = n_nws_pts;
            for (p = 0; p < n_local_nws_pts; p++) {
                nws_pts(n_nws_pts++) = local_nws_pts(p);
            }
            for (q = 0; q < n_local_nws_pts - 1; q++) {
                nws_seg(0, n_nws_seg) = map + q;
                nws_seg(1, n_nws_seg) = map + q + 1;
                n_nws_seg++;
            }

            /////////////////////////////////////////
            ////// CONSTRUCT THE nw SURFACE /////////
            /////////////////////////////////////////
            if (n_local_nws_pts > 0) {
                EDGE(F, vF, S, i, j, k, Nx, Ny, Nz, nw_pts, n_nw_pts, nw_tris,
                     n_nw_tris, local_nws_pts, n_local_nws_pts);
            } else {
                MC(F, vF, S, i, j, k, nw_pts, n_nw_pts, nw_tris, n_nw_tris);
            }
        }

        /////////////////////////////////////////
        ////// CONSTRUCT THE nw SURFACE /////////
        /////////////////////////////////////////

        else if (Fluid_Interface(F, S, vF, i, j, k) == 1) {
            MC(F, vF, S, i, j, k, nw_pts, n_nw_pts, nw_tris, n_nw_tris);
        }
        //******END OF BLOB PMMC*********************************************

        //*******************************************************************
        // Compute the Interfacial Areas, Common Line length for blob p
        // nw surface
        for (r = n_nw_tris_beg; r < n_nw_tris; r++) {
            A = nw_pts(nw_tris(0, r));
            B = nw_pts(nw_tris(1, r));
            C = nw_pts(nw_tris(2, r));
            // Compute length of sides (assume dx=dy=dz)
            s1 = sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y) +
                      (A.z - B.z) * (A.z - B.z));
            s2 = sqrt((A.x - C.x) * (A.x - C.x) + (A.y - C.y) * (A.y - C.y) +
                      (A.z - C.z) * (A.z - C.z));
            s3 = sqrt((B.x - C.x) * (B.x - C.x) + (B.y - C.y) * (B.y - C.y) +
                      (B.z - C.z) * (B.z - C.z));
            s = 0.5 * (s1 + s2 + s3);
            awn = awn + sqrt(s * (s - s1) * (s - s2) * (s - s3));
        }
        for (r = n_ns_tris_beg; r < n_ns_tris; r++) {
            A = ns_pts(ns_tris(0, r));
            B = ns_pts(ns_tris(1, r));
            C = ns_pts(ns_tris(2, r));
            // Compute length of sides (assume dx=dy=dz)
            s1 = sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y) +
                      (A.z - B.z) * (A.z - B.z));
            s2 = sqrt((A.x - C.x) * (A.x - C.x) + (A.y - C.y) * (A.y - C.y) +
                      (A.z - C.z) * (A.z - C.z));
            s3 = sqrt((B.x - C.x) * (B.x - C.x) + (B.y - C.y) * (B.y - C.y) +
                      (B.z - C.z) * (B.z - C.z));
            s = 0.5 * (s1 + s2 + s3);
            ans = ans + sqrt(s * (s - s1) * (s - s2) * (s - s3));
        }
        for (r = n_ws_tris_beg; r < n_ws_tris; r++) {
            A = ws_pts(ws_tris(0, r));
            B = ws_pts(ws_tris(1, r));
            C = ws_pts(ws_tris(2, r));
            // Compute length of sides (assume dx=dy=dz)
            s1 = sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y) +
                      (A.z - B.z) * (A.z - B.z));
            s2 = sqrt((A.x - C.x) * (A.x - C.x) + (A.y - C.y) * (A.y - C.y) +
                      (A.z - C.z) * (A.z - C.z));
            s3 = sqrt((B.x - C.x) * (B.x - C.x) + (B.y - C.y) * (B.y - C.y) +
                      (B.z - C.z) * (B.z - C.z));
            s = 0.5 * (s1 + s2 + s3);
            aws = aws + sqrt(s * (s - s1) * (s - s2) * (s - s3));
        }
        //*******************************************************************
        // Reset the triangle counts to zero
        n_nw_pts = 0, n_ns_pts = 0, n_ws_pts = 0, n_nws_pts = 0, map = 0;
        n_nw_tris = 0, n_ns_tris = 0, n_ws_tris = 0, n_nws_seg = 0;

        n_nw_tris_beg = n_nw_tris;
        n_ns_tris_beg = n_ns_tris;
        n_ws_tris_beg = n_ws_tris;
        //		n_nws_seg_beg = n_nws_seg;
        //*******************************************************************
    }
}
//--------------------------------------------------------------------------------------------------------
inline void pmmc_ConstructLocalCube(
    DoubleArray &SignDist, DoubleArray &Phase, double solid_isovalue,
    double fluid_isovalue, DTMutableList<Point> &nw_pts, IntArray &nw_tris,
    DoubleArray &values, DTMutableList<Point> &ns_pts, IntArray &ns_tris,
    DTMutableList<Point> &ws_pts, IntArray &ws_tris,
    DTMutableList<Point> &local_nws_pts, DTMutableList<Point> &nws_pts,
    IntArray &nws_seg, DTMutableList<Point> &local_sol_pts,
    IntArray &local_sol_tris, int &n_local_sol_tris, int &n_local_sol_pts,
    int &n_nw_pts, int &n_nw_tris, int &n_ws_pts, int &n_ws_tris,
    int &n_ns_tris, int &n_ns_pts, int &n_local_nws_pts, int &n_nws_pts,
    int &n_nws_seg, int i, int j, int k, int Nx, int Ny, int Nz) {

    int p, q, map;
    Point A, B, C, P;

    // Only the local values are constructed and retained! (set counts to zero to force this)
    n_local_sol_tris = 0;
    n_local_sol_pts = 0;
    n_local_nws_pts = 0;

    n_nw_pts = 0, n_ns_pts = 0, n_ws_pts = 0, n_nws_pts = 0, map = 0;
    n_nw_tris = 0, n_ns_tris = 0, n_ws_tris = 0, n_nws_seg = 0;

    // if there is a solid phase interface in the grid cell
    if (Interface(SignDist, solid_isovalue, i, j, k) == 1) {

        /////////////////////////////////////////
        /// CONSTRUCT THE LOCAL SOLID SURFACE ///
        /////////////////////////////////////////

        // find the local solid surface using the regular Marching Cubes algorithm
        SolidMarchingCubes(SignDist, 0.0, Phase, fluid_isovalue, i, j, k, Nx,
                           Ny, Nz, local_sol_pts, n_local_sol_pts,
                           local_sol_tris, n_local_sol_tris, values);

        /////////////////////////////////////////
        //////// TRIM THE SOLID SURFACE /////////
        /////////////////////////////////////////
        TRIM(local_sol_pts, n_local_sol_pts, fluid_isovalue, local_sol_tris,
             n_local_sol_tris, ns_pts, n_ns_pts, ns_tris, n_ns_tris, ws_pts,
             n_ws_pts, ws_tris, n_ws_tris, values, local_nws_pts,
             n_local_nws_pts);

        /////////////////////////////////////////
        //////// WRITE COMMON LINE POINTS ///////
        ////////      TO MAIN ARRAYS      ///////
        /////////////////////////////////////////
        // SORT THE LOCAL COMMON LINE POINTS
        /////////////////////////////////////////
        // Make sure the first common line point is on a face
        // Common curve points are located pairwise and must
        // be searched and rearranged accordingly
        if (n_local_nws_pts > 0) {
            for (p = 0; p < n_local_nws_pts - 1; p++) {
                P = local_nws_pts(p);
                if (P.x == 1.0 * i || P.x == 1.0 * (i + 1) || P.y == 1.0 * j ||
                    P.y == 1.0 * (j + 1) || P.z == 1.0 * k ||
                    P.z == 1.0 * (k + 1)) {
                    if (p % 2 == 0) {
                        // even points
                        // Swap the pair of points
                        local_nws_pts(p) = local_nws_pts(0);
                        local_nws_pts(0) = P;
                        P = local_nws_pts(p + 1);
                        local_nws_pts(p + 1) = local_nws_pts(1);
                        local_nws_pts(1) = P;
                        p = n_local_nws_pts;

                    } else {
                        // odd points - flip the order
                        local_nws_pts(p) = local_nws_pts(p - 1);
                        local_nws_pts(p - 1) = P;
                        p -= 2;
                    }
                    // guarantee exit from the loop
                }
            }
            // Two common curve points per triangle
            // 0-(1=2)-(3=4)-...
            for (p = 1; p < n_local_nws_pts - 1; p += 2) {
                A = local_nws_pts(p);
                for (q = p + 1; q < n_local_nws_pts; q++) {
                    B = local_nws_pts(q);
                    if (A.x == B.x && A.y == B.y && A.z == B.z) {
                        if (q % 2 == 0) {
                            // even points
                            // Swap the pair of points
                            local_nws_pts(q) = local_nws_pts(p + 1);
                            local_nws_pts(p + 1) = B;
                            B = local_nws_pts(q + 1);
                            local_nws_pts(q + 1) = local_nws_pts(p + 2);
                            local_nws_pts(p + 2) = B;
                            q = n_local_nws_pts;

                        } else {
                            // odd points - flip the order
                            local_nws_pts(q) = local_nws_pts(q - 1);
                            local_nws_pts(q - 1) = B;
                            q -= 2;
                        }
                    }
                }
            }
            map = n_nws_pts = 0;
            nws_pts(n_nws_pts++) = local_nws_pts(0);
            for (p = 2; p < n_local_nws_pts; p += 2) {
                nws_pts(n_nws_pts++) = local_nws_pts(p);
            }
            nws_pts(n_nws_pts++) = local_nws_pts(n_local_nws_pts - 1);

            for (q = 0; q < n_nws_pts - 1; q++) {
                nws_seg(0, n_nws_seg) = map + q;
                nws_seg(1, n_nws_seg) = map + q + 1;
                n_nws_seg++;
            }
            // End of the common line sorting algorithm
            /////////////////////////////////////////
        }

        /////////////////////////////////////////
        ////// CONSTRUCT THE nw SURFACE /////////
        /////////////////////////////////////////
        if (n_local_nws_pts > 0) {
            n_nw_tris = 0;
            EDGE(Phase, fluid_isovalue, SignDist, i, j, k, Nx, Ny, Nz, nw_pts,
                 n_nw_pts, nw_tris, n_nw_tris, nws_pts, n_nws_pts);
        } else {
            MC(Phase, fluid_isovalue, SignDist, i, j, k, nw_pts, n_nw_pts,
               nw_tris, n_nw_tris);
        }
    }

    /////////////////////////////////////////
    ////// CONSTRUCT THE nw SURFACE /////////
    /////////////////////////////////////////

    else if (Fluid_Interface(Phase, SignDist, fluid_isovalue, i, j, k) == 1) {
        MC(Phase, fluid_isovalue, SignDist, i, j, k, nw_pts, n_nw_pts, nw_tris,
           n_nw_tris);
    }
}
//--------------------------------------------------------------------------------------------------------
inline void pmmc_MeshGradient(DoubleArray &f, DoubleArray &fx, DoubleArray &fy,
                              DoubleArray &fz, int Nx, int Ny, int Nz) {
    int i, j, k;
    // Compute the Gradient everywhere except the halo region
    for (k = 1; k < Nz - 1; k++) {
        for (j = 1; j < Ny - 1; j++) {
            for (i = 1; i < Nx - 1; i++) {
                // Compute all of the derivatives using finite differences
                fx(i, j, k) = 0.5 * (f(i + 1, j, k) - f(i - 1, j, k));
                fy(i, j, k) = 0.5 * (f(i, j + 1, k) - f(i, j - 1, k));
                fz(i, j, k) = 0.5 * (f(i, j, k + 1) - f(i, j, k - 1));
            }
        }
    }
}
//--------------------------------------------------------------------------------------------------------
inline void pmmc_MeshCurvature(DoubleArray &f, DoubleArray &MeanCurvature,
                               DoubleArray &GaussCurvature, int Nx, int Ny,
                               int Nz) {
    // Mesh spacing is taken to be one to simplify the calculation
    int i, j, k;
    double denominator;
    double fxx, fyy, fzz, fxy, fxz, fyz, fx, fy, fz;

    // Compute the curvature everywhere except the halo region
    for (k = 1; k < Nz - 1; k++) {
        for (j = 1; j < Ny - 1; j++) {
            for (i = 1; i < Nx - 1; i++) {
                // Compute all of the derivatives using finite differences
                fx = 0.5 * (f(i + 1, j, k) - f(i - 1, j, k));
                fy = 0.5 * (f(i, j + 1, k) - f(i, j - 1, k));
                fz = 0.5 * (f(i, j, k + 1) - f(i, j, k - 1));
                fxx = f(i + 1, j, k) - 2.0 * f(i, j, k) + f(i - 1, j, k);
                fyy = f(i, j + 1, k) - 2.0 * f(i, j, k) + f(i, j - 1, k);
                fzz = f(i, j, k + 1) - 2.0 * f(i, j, k) + f(i, j, k - 1);
                fxy = 0.25 * (f(i + 1, j + 1, k) - f(i + 1, j - 1, k) -
                              f(i - 1, j + 1, k) + f(i - 1, j - 1, k));
                fxz = 0.25 * (f(i + 1, j, k + 1) - f(i + 1, j, k - 1) -
                              f(i - 1, j, k + 1) + f(i - 1, j, k - 1));
                fyz = 0.25 * (f(i, j + 1, k + 1) - f(i, j + 1, k - 1) -
                              f(i, j - 1, k + 1) + f(i, j - 1, k - 1));
                // Evaluate the Mean Curvature
                denominator = pow(sqrt(fx * fx + fy * fy + fz * fz), 3);
                if (denominator == 0.0) {
                    MeanCurvature(i, j, k) = 0.0;
                } else {
                    MeanCurvature(i, j, k) =
                        (1.0 / denominator) *
                        ((fyy + fzz) * fx * fx + (fxx + fzz) * fy * fy +
                         (fxx + fyy) * fz * fz - 2.0 * fx * fy * fxy -
                         2.0 * fx * fz * fxz - 2.0 * fy * fz * fyz);
                }
                // Evaluate the Gaussian Curvature
                denominator = pow(fx * fx + fy * fy + fz * fz, 2);
                if (denominator == 0.0) {
                    GaussCurvature(i, j, k) = 0.0;
                } else {

                    GaussCurvature(i, j, k) =
                        (1.0 / denominator) *
                        (fx * fx * (fyy * fzz - fyz * fyz) +
                         fy * fy * (fxx * fzz - fxz * fxz) +
                         fz * fz * (fxx * fyy - fxy * fxy) +
                         2.0 * (fx * fy * (fxz * fyz - fxy * fzz) +
                                fy * fz * (fxy * fxz - fyz * fxx) +
                                fx * fz * (fxy * fyz - fxz * fyy)));
                }
            }
        }
    }
}
//--------------------------------------------------------------------------------------------------------
inline int pmmc_CubeListFromMesh(IntArray &cubeList, int ncubes, int Nx, int Ny,
                                 int Nz) {
    NULL_USE(ncubes);

    int i, j, k, nc;
    nc = 0;
    //...........................................................................
    // Set up the cube list (very regular in this case due to lack of blob-ID)
    for (k = 1; k < Nz - 2; k++) {
        for (j = 1; j < Ny - 2; j++) {
            for (i = 1; i < Nx - 2; i++) {
                cubeList(0, nc) = i;
                cubeList(1, nc) = j;
                cubeList(2, nc) = k;
                nc++;
            }
        }
    }
    return nc;
}
//--------------------------------------------------------------------------------------------------------
inline void pmmc_CubeListFromBlobs() {}
//--------------------------------------------------------------------------------------------------------
inline double pmmc_CubeSurfaceArea(DTMutableList<Point> &Points,
                                   IntArray &Triangles, int ntris) {
    int r;
    double temp, area, s, s1, s2, s3;
    Point A, B, C;
    area = 0.0;
    for (r = 0; r < ntris; r++) {
        A = Points(Triangles(0, r));
        B = Points(Triangles(1, r));
        C = Points(Triangles(2, r));
        // Compute length of sides (assume dx=dy=dz)
        s1 = sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y) +
                  (A.z - B.z) * (A.z - B.z));
        s2 = sqrt((A.x - C.x) * (A.x - C.x) + (A.y - C.y) * (A.y - C.y) +
                  (A.z - C.z) * (A.z - C.z));
        s3 = sqrt((B.x - C.x) * (B.x - C.x) + (B.y - C.y) * (B.y - C.y) +
                  (B.z - C.z) * (B.z - C.z));
        s = 0.5 * (s1 + s2 + s3);
        temp = s * (s - s1) * (s - s2) * (s - s3);
        if (temp > 0.0)
            area += sqrt(temp);
    }
    return area;
}
//--------------------------------------------------------------------------------------------------------
inline double pmmc_CubeCurveLength(DTMutableList<Point> &Points, int npts) {
    int p;
    double s, lwns;
    Point A, B;
    lwns = 0.0;
    for (p = 0; p < npts - 1; p++) {
        // Extract the line segment
        A = Points(p);
        B = Points(p + 1);
        // Compute the length of the segment
        s = sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y) +
                 (A.z - B.z) * (A.z - B.z));
        // Add the length to the common line
        lwns += s;
    }
    return lwns;
}
//--------------------------------------------------------------------------------------------------------
inline void pmmc_CubeTrimSurfaceInterpValues(
    DoubleArray &CubeValues, DoubleArray &MeshValues, DoubleArray &SignDist,
    DTMutableList<Point> &Points, IntArray &Triangles,
    DoubleArray &SurfaceValues, DoubleArray &DistanceValues, int i, int j,
    int k, int npts, int ntris, double mindist, double &area,
    double &integral) {
    // mindist - minimum distance to consider in the average (in voxel lengths)
    Point A, B, C;
    int p;
    double vA, vB, vC;
    double dA, dB, dC;
    double x, y, z;
    double s, s1, s2, s3, temp;
    double a, b, c, d, e, f, g, h;

    // Copy the curvature values for the cube
    CubeValues(0, 0, 0) = MeshValues(i, j, k);
    CubeValues(1, 0, 0) = MeshValues(i + 1, j, k);
    CubeValues(0, 1, 0) = MeshValues(i, j + 1, k);
    CubeValues(1, 1, 0) = MeshValues(i + 1, j + 1, k);
    CubeValues(0, 0, 1) = MeshValues(i, j, k + 1);
    CubeValues(1, 0, 1) = MeshValues(i + 1, j, k + 1);
    CubeValues(0, 1, 1) = MeshValues(i, j + 1, k + 1);
    CubeValues(1, 1, 1) = MeshValues(i + 1, j + 1, k + 1);

    // trilinear coefficients: f(x,y,z) = a+bx+cy+dz+exy+fxz+gyz+hxyz
    // Evaluate the coefficients
    a = CubeValues(0, 0, 0);
    b = CubeValues(1, 0, 0) - a;
    c = CubeValues(0, 1, 0) - a;
    d = CubeValues(0, 0, 1) - a;
    e = CubeValues(1, 1, 0) - a - b - c;
    f = CubeValues(1, 0, 1) - a - b - d;
    g = CubeValues(0, 1, 1) - a - c - d;
    h = CubeValues(1, 1, 1) - a - b - c - d - e - f - g;

    for (p = 0; p < npts; p++) {
        A = Points(p);
        x = A.x - 1.0 * i;
        y = A.y - 1.0 * j;
        z = A.z - 1.0 * k;
        SurfaceValues(p) = a + b * x + c * y + d * z + e * x * y + f * x * z +
                           g * y * z + h * x * y * z;
    }

    // Copy the signed distance values for the cube
    CubeValues(0, 0, 0) = SignDist(i, j, k);
    CubeValues(1, 0, 0) = SignDist(i + 1, j, k);
    CubeValues(0, 1, 0) = SignDist(i, j + 1, k);
    CubeValues(1, 1, 0) = SignDist(i + 1, j + 1, k);
    CubeValues(0, 0, 1) = SignDist(i, j, k + 1);
    CubeValues(1, 0, 1) = SignDist(i + 1, j, k + 1);
    CubeValues(0, 1, 1) = SignDist(i, j + 1, k + 1);
    CubeValues(1, 1, 1) = SignDist(i + 1, j + 1, k + 1);

    // trilinear coefficients: f(x,y,z) = a+bx+cy+dz+exy+fxz+gyz+hxyz
    // Evaluate the coefficients
    a = CubeValues(0, 0, 0);
    b = CubeValues(1, 0, 0) - a;
    c = CubeValues(0, 1, 0) - a;
    d = CubeValues(0, 0, 1) - a;
    e = CubeValues(1, 1, 0) - a - b - c;
    f = CubeValues(1, 0, 1) - a - b - d;
    g = CubeValues(0, 1, 1) - a - c - d;
    h = CubeValues(1, 1, 1) - a - b - c - d - e - f - g;

    for (p = 0; p < npts; p++) {
        A = Points(p);
        x = A.x - 1.0 * i;
        y = A.y - 1.0 * j;
        z = A.z - 1.0 * k;
        DistanceValues(p) = a + b * x + c * y + d * z + e * x * y + f * x * z +
                            g * y * z + h * x * y * z;
    }

    for (int r = 0; r < ntris; r++) {
        A = Points(Triangles(0, r));
        B = Points(Triangles(1, r));
        C = Points(Triangles(2, r));

        vA = SurfaceValues(Triangles(0, r));
        vB = SurfaceValues(Triangles(1, r));
        vC = SurfaceValues(Triangles(2, r));

        dA = DistanceValues(Triangles(0, r));
        dB = DistanceValues(Triangles(1, r));
        dC = DistanceValues(Triangles(2, r));

        // Compute length of sides (assume dx=dy=dz)
        s1 = sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y) +
                  (A.z - B.z) * (A.z - B.z));
        s2 = sqrt((A.x - C.x) * (A.x - C.x) + (A.y - C.y) * (A.y - C.y) +
                  (A.z - C.z) * (A.z - C.z));
        s3 = sqrt((B.x - C.x) * (B.x - C.x) + (B.y - C.y) * (B.y - C.y) +
                  (B.z - C.z) * (B.z - C.z));
        s = 0.5 * (s1 + s2 + s3);
        temp = s * (s - s1) * (s - s2) * (s - s3);

        if (temp > 0.0) {
            if (dA > mindist) {
                integral += sqrt(temp) * 0.33333333333333333 * (vA);
                area += sqrt(temp) * 0.33333333333333333;
            }
            if (dB > mindist) {
                integral += sqrt(temp) * 0.33333333333333333 * (vB);
                area += sqrt(temp) * 0.33333333333333333;
            }
            if (dC > mindist) {
                integral += sqrt(temp) * 0.33333333333333333 * (vC);
                area += sqrt(temp) * 0.33333333333333333;
            }
        }
    }
}
inline void pmmc_CubeTrimSurfaceInterpInverseValues(
    DoubleArray &CubeValues, DoubleArray &MeshValues, DoubleArray &SignDist,
    DTMutableList<Point> &Points, IntArray &Triangles,
    DoubleArray &SurfaceValues, DoubleArray &DistanceValues, int i, int j,
    int k, int npts, int ntris, double mindist, double &area,
    double &integral) {
    // mindist - minimum distance to consider in the average (in voxel lengths)
    Point A, B, C;
    int p;
    double vA, vB, vC;
    double dA, dB, dC;
    double x, y, z;
    double s, s1, s2, s3, temp;
    double a, b, c, d, e, f, g, h;

    // Copy the curvature values for the cube
    CubeValues(0, 0, 0) = MeshValues(i, j, k);
    CubeValues(1, 0, 0) = MeshValues(i + 1, j, k);
    CubeValues(0, 1, 0) = MeshValues(i, j + 1, k);
    CubeValues(1, 1, 0) = MeshValues(i + 1, j + 1, k);
    CubeValues(0, 0, 1) = MeshValues(i, j, k + 1);
    CubeValues(1, 0, 1) = MeshValues(i + 1, j, k + 1);
    CubeValues(0, 1, 1) = MeshValues(i, j + 1, k + 1);
    CubeValues(1, 1, 1) = MeshValues(i + 1, j + 1, k + 1);

    // trilinear coefficients: f(x,y,z) = a+bx+cy+dz+exy+fxz+gyz+hxyz
    // Evaluate the coefficients
    a = CubeValues(0, 0, 0);
    b = CubeValues(1, 0, 0) - a;
    c = CubeValues(0, 1, 0) - a;
    d = CubeValues(0, 0, 1) - a;
    e = CubeValues(1, 1, 0) - a - b - c;
    f = CubeValues(1, 0, 1) - a - b - d;
    g = CubeValues(0, 1, 1) - a - c - d;
    h = CubeValues(1, 1, 1) - a - b - c - d - e - f - g;

    for (p = 0; p < npts; p++) {
        A = Points(p);
        x = A.x - 1.0 * i;
        y = A.y - 1.0 * j;
        z = A.z - 1.0 * k;
        SurfaceValues(p) = a + b * x + c * y + d * z + e * x * y + f * x * z +
                           g * y * z + h * x * y * z;
    }

    // Copy the signed distance values for the cube
    CubeValues(0, 0, 0) = SignDist(i, j, k);
    CubeValues(1, 0, 0) = SignDist(i + 1, j, k);
    CubeValues(0, 1, 0) = SignDist(i, j + 1, k);
    CubeValues(1, 1, 0) = SignDist(i + 1, j + 1, k);
    CubeValues(0, 0, 1) = SignDist(i, j, k + 1);
    CubeValues(1, 0, 1) = SignDist(i + 1, j, k + 1);
    CubeValues(0, 1, 1) = SignDist(i, j + 1, k + 1);
    CubeValues(1, 1, 1) = SignDist(i + 1, j + 1, k + 1);

    // trilinear coefficients: f(x,y,z) = a+bx+cy+dz+exy+fxz+gyz+hxyz
    // Evaluate the coefficients
    a = CubeValues(0, 0, 0);
    b = CubeValues(1, 0, 0) - a;
    c = CubeValues(0, 1, 0) - a;
    d = CubeValues(0, 0, 1) - a;
    e = CubeValues(1, 1, 0) - a - b - c;
    f = CubeValues(1, 0, 1) - a - b - d;
    g = CubeValues(0, 1, 1) - a - c - d;
    h = CubeValues(1, 1, 1) - a - b - c - d - e - f - g;

    for (p = 0; p < npts; p++) {
        A = Points(p);
        x = A.x - 1.0 * i;
        y = A.y - 1.0 * j;
        z = A.z - 1.0 * k;
        DistanceValues(p) = a + b * x + c * y + d * z + e * x * y + f * x * z +
                            g * y * z + h * x * y * z;
    }

    for (int r = 0; r < ntris; r++) {
        A = Points(Triangles(0, r));
        B = Points(Triangles(1, r));
        C = Points(Triangles(2, r));

        vA = SurfaceValues(Triangles(0, r));
        vB = SurfaceValues(Triangles(1, r));
        vC = SurfaceValues(Triangles(2, r));

        dA = DistanceValues(Triangles(0, r));
        dB = DistanceValues(Triangles(1, r));
        dC = DistanceValues(Triangles(2, r));

        // Compute length of sides (assume dx=dy=dz)
        s1 = sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y) +
                  (A.z - B.z) * (A.z - B.z));
        s2 = sqrt((A.x - C.x) * (A.x - C.x) + (A.y - C.y) * (A.y - C.y) +
                  (A.z - C.z) * (A.z - C.z));
        s3 = sqrt((B.x - C.x) * (B.x - C.x) + (B.y - C.y) * (B.y - C.y) +
                  (B.z - C.z) * (B.z - C.z));
        s = 0.5 * (s1 + s2 + s3);
        temp = s * (s - s1) * (s - s2) * (s - s3);

        if (temp > 0.0) {
            if (dA > mindist && vA != 0.0) {
                integral += sqrt(temp) * 0.33333333333333333 * (1.0 / vA);
                area += sqrt(temp) * 0.33333333333333333;
            }
            if (dB > mindist && vB != 0.0) {
                integral += sqrt(temp) * 0.33333333333333333 * (1.0 / vB);
                area += sqrt(temp) * 0.33333333333333333;
            }
            if (dC > mindist && vC != 0.0) {
                integral += sqrt(temp) * 0.33333333333333333 * (1.0 / vC);
                area += sqrt(temp) * 0.33333333333333333;
            }
        }
    }
}

inline double pmmc_CubeSurfaceInterpValue(DoubleArray &CubeValues,
                                          DoubleArray &MeshValues,
                                          DTMutableList<Point> &Points,
                                          IntArray &Triangles,
                                          DoubleArray &SurfaceValues, int i,
                                          int j, int k, int npts, int ntris) {
    Point A, B, C;
    int p;
    double vA, vB, vC;
    double x, y, z;
    double s, s1, s2, s3, temp;
    double a, b, c, d, e, f, g, h;
    double integral;

    // Copy the curvature values for the cube
    CubeValues(0, 0, 0) = MeshValues(i, j, k);
    CubeValues(1, 0, 0) = MeshValues(i + 1, j, k);
    CubeValues(0, 1, 0) = MeshValues(i, j + 1, k);
    CubeValues(1, 1, 0) = MeshValues(i + 1, j + 1, k);
    CubeValues(0, 0, 1) = MeshValues(i, j, k + 1);
    CubeValues(1, 0, 1) = MeshValues(i + 1, j, k + 1);
    CubeValues(0, 1, 1) = MeshValues(i, j + 1, k + 1);
    CubeValues(1, 1, 1) = MeshValues(i + 1, j + 1, k + 1);

    // trilinear coefficients: f(x,y,z) = a+bx+cy+dz+exy+fxz+gyz+hxyz
    // Evaluate the coefficients
    a = CubeValues(0, 0, 0);
    b = CubeValues(1, 0, 0) - a;
    c = CubeValues(0, 1, 0) - a;
    d = CubeValues(0, 0, 1) - a;
    e = CubeValues(1, 1, 0) - a - b - c;
    f = CubeValues(1, 0, 1) - a - b - d;
    g = CubeValues(0, 1, 1) - a - c - d;
    h = CubeValues(1, 1, 1) - a - b - c - d - e - f - g;

    for (p = 0; p < npts; p++) {
        A = Points(p);
        x = A.x - 1.0 * i;
        y = A.y - 1.0 * j;
        z = A.z - 1.0 * k;
        SurfaceValues(p) = a + b * x + c * y + d * z + e * x * y + f * x * z +
                           g * y * z + h * x * y * z;
    }

    integral = 0.0;
    for (int r = 0; r < ntris; r++) {
        A = Points(Triangles(0, r));
        B = Points(Triangles(1, r));
        C = Points(Triangles(2, r));
        vA = SurfaceValues(Triangles(0, r));
        vB = SurfaceValues(Triangles(1, r));
        vC = SurfaceValues(Triangles(2, r));
        // Compute length of sides (assume dx=dy=dz)
        s1 = sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y) +
                  (A.z - B.z) * (A.z - B.z));
        s2 = sqrt((A.x - C.x) * (A.x - C.x) + (A.y - C.y) * (A.y - C.y) +
                  (A.z - C.z) * (A.z - C.z));
        s3 = sqrt((B.x - C.x) * (B.x - C.x) + (B.y - C.y) * (B.y - C.y) +
                  (B.z - C.z) * (B.z - C.z));
        s = 0.5 * (s1 + s2 + s3);
        temp = s * (s - s1) * (s - s2) * (s - s3);
        if (temp > 0.0)
            integral += sqrt(temp) * 0.33333333333333333 * (vA + vB + vC);
    }
    return integral;
}
//--------------------------------------------------------------------------------------------------------
inline double pmmc_CubeCurveInterpValue(DoubleArray &CubeValues,
                                        DoubleArray &CurveValues,
                                        DTMutableList<Point> &Points, int i,
                                        int j, int k, int npts) {
    int p;
    Point A, B;
    double vA, vB;
    double x, y, z;
    //	double s,s1,s2,s3,temp;
    double a, b, c, d, e, f, g, h;
    double integral;
    double length;

    // trilinear coefficients: f(x,y,z) = a+bx+cy+dz+exy+fxz+gyz+hxyz
    // Evaluate the coefficients
    a = CubeValues(0, 0, 0);
    b = CubeValues(1, 0, 0) - a;
    c = CubeValues(0, 1, 0) - a;
    d = CubeValues(0, 0, 1) - a;
    e = CubeValues(1, 1, 0) - a - b - c;
    f = CubeValues(1, 0, 1) - a - b - d;
    g = CubeValues(0, 1, 1) - a - c - d;
    h = CubeValues(1, 1, 1) - a - b - c - d - e - f - g;

    for (p = 0; p < npts; p++) {
        A = Points(p);
        x = A.x - 1.0 * i;
        y = A.y - 1.0 * j;
        z = A.z - 1.0 * k;
        CurveValues(p) = a + b * x + c * y + d * z + e * x * y + f * x * z +
                         g * y * z + h * x * y * z;
    }

    integral = 0.0;
    for (p = 0; p < npts - 1; p++) {
        // Extract the line segment
        A = Points(p);
        B = Points(p + 1);
        vA = CurveValues(p);
        vB = CurveValues(p + 1);
        // Compute the length of the segment
        length = sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y) +
                      (A.z - B.z) * (A.z - B.z));
        integral += 0.5 * length * (vA + vB);
    }
    return integral;
}
//--------------------------------------------------------------------------------------------------------
inline double pmmc_CubeContactAngle(DoubleArray &CubeValues,
                                    DoubleArray &CurveValues, DoubleArray &Fx,
                                    DoubleArray &Fy, DoubleArray &Fz,
                                    DoubleArray &Sx, DoubleArray &Sy,
                                    DoubleArray &Sz,
                                    DTMutableList<Point> &Points, int i, int j,
                                    int k, int npts) {
    int p;
    Point A, B;
    double vA, vB;
    double x, y, z;
    //	double s,s1,s2,s3,temp;
    double a, b, c, d, e, f, g, h;
    double integral;
    double length;
    double denom;

    // theta = acos ( -(gradF*gradS) / (|gradF| |gradS|) )

    denom =
        sqrt(pow(Fx(i, j, k), 2) + pow(Fy(i, j, k), 2) + pow(Fz(i, j, k), 2)) *
        sqrt(pow(Sx(i, j, k), 2) + pow(Sy(i, j, k), 2) + pow(Sz(i, j, k), 2));
    if (denom == 0.0)
        denom = 1.0;
    CubeValues(0, 0, 0) =
        -(Fx(i, j, k) * Sx(i, j, k) + Fy(i, j, k) * Sy(i, j, k) +
          Fz(i, j, k) * Sz(i, j, k)) /
        (denom);

    denom = sqrt(pow(Fx(i + 1, j, k), 2) + pow(Fy(i + 1, j, k), 2) +
                 pow(Fz(i + 1, j, k), 2)) *
            sqrt(pow(Sx(i + 1, j, k), 2) + pow(Sy(i + 1, j, k), 2) +
                 pow(Sz(i + 1, j, k), 2));
    if (denom == 0.0)
        denom = 1.0;
    CubeValues(1, 0, 0) = -(Fx(i + 1, j, k) * Sx(i + 1, j, k) +
                            Fy(i + 1, j, k) * Sy(i + 1, j, k) +
                            Fz(i + 1, j, k) * Sz(i + 1, j, k)) /
                          (denom);

    denom = sqrt(pow(Fx(i, j + 1, k), 2) + pow(Fy(i, j + 1, k), 2) +
                 pow(Fz(i, j + 1, k), 2)) *
            sqrt(pow(Sx(i, j + 1, k), 2) + pow(Sy(i, j + 1, k), 2) +
                 pow(Sz(i, j + 1, k), 2));
    if (denom == 0.0)
        denom = 1.0;
    CubeValues(0, 1, 0) = -(Fx(i, j + 1, k) * Sx(i, j + 1, k) +
                            Fy(i, j + 1, k) * Sy(i, j + 1, k) +
                            Fz(i, j + 1, k) * Sz(i, j + 1, k)) /
                          (denom);

    denom = sqrt(pow(Fx(i, j, k + 1), 2) + pow(Fy(i, j, k + 1), 2) +
                 pow(Fz(i, j, k + 1), 2)) *
            sqrt(pow(Sx(i, j, k + 1), 2) + pow(Sy(i, j, k + 1), 2) +
                 pow(Sz(i, j, k + 1), 2));
    if (denom == 0.0)
        denom = 1.0;
    CubeValues(0, 0, 1) = -(Fx(i, j, k + 1) * Sx(i, j, k + 1) +
                            Fy(i, j, k + 1) * Sy(i, j, k + 1) +
                            Fz(i, j, k + 1) * Sz(i, j, k + 1)) /
                          (denom);

    denom = sqrt(pow(Fx(i + 1, j + 1, k), 2) + pow(Fy(i + 1, j + 1, k), 2) +
                 pow(Fz(i + 1, j + 1, k), 2)) *
            sqrt(pow(Sx(i + 1, j + 1, k), 2) + pow(Sy(i + 1, j + 1, k), 2) +
                 pow(Sz(i + 1, j + 1, k), 2));
    if (denom == 0.0)
        denom = 1.0;
    CubeValues(1, 1, 0) = -(Fx(i + 1, j + 1, k) * Sx(i + 1, j + 1, k) +
                            Fy(i + 1, j + 1, k) * Sy(i + 1, j + 1, k) +
                            Fz(i + 1, j + 1, k) * Sz(i + 1, j + 1, k)) /
                          (denom);

    denom = sqrt(pow(Fx(i + 1, j, k + 1), 2) + pow(Fy(i + 1, j, k + 1), 2) +
                 pow(Fz(i + 1, j, k + 1), 2)) *
            sqrt(pow(Sx(i + 1, j, k + 1), 2) + pow(Sy(i + 1, j, k + 1), 2) +
                 pow(Sz(i + 1, j, k + 1), 2));
    if (denom == 0.0)
        denom = 1.0;
    CubeValues(1, 0, 1) = -(Fx(i + 1, j, k + 1) * Sx(i + 1, j, k + 1) +
                            Fy(i + 1, j, k + 1) * Sy(i + 1, j, k + 1) +
                            Fz(i + 1, j, k + 1) * Sz(i + 1, j, k + 1)) /
                          (denom);

    denom = sqrt(pow(Fx(i, j + 1, k + 1), 2) + pow(Fy(i, j + 1, k + 1), 2) +
                 pow(Fz(i, j + 1, k + 1), 2)) *
            sqrt(pow(Sx(i, j + 1, k + 1), 2) + pow(Sy(i, j + 1, k + 1), 2) +
                 pow(Sz(i, j + 1, k + 1), 2));
    if (denom == 0.0)
        denom = 1.0;
    CubeValues(0, 1, 1) = -(Fx(i, j + 1, k + 1) * Sx(i, j + 1, k + 1) +
                            Fy(i, j + 1, k + 1) * Sy(i, j + 1, k + 1) +
                            Fz(i, j + 1, k + 1) * Sz(i, j + 1, k + 1)) /
                          (denom);

    denom =
        sqrt(pow(Fx(i + 1, j + 1, k + 1), 2) + pow(Fy(i + 1, j + 1, k + 1), 2) +
             pow(Fz(i + 1, j + 1, k + 1), 2)) *
        sqrt(pow(Sx(i + 1, j + 1, k + 1), 2) + pow(Sy(i + 1, j + 1, k + 1), 2) +
             pow(Sz(i + 1, j + 1, k + 1), 2));
    if (denom == 0.0)
        denom = 1.0;
    CubeValues(1, 1, 1) = -(Fx(i + 1, j + 1, k + 1) * Sx(i + 1, j + 1, k + 1) +
                            Fy(i + 1, j + 1, k + 1) * Sy(i + 1, j + 1, k + 1) +
                            Fz(i + 1, j + 1, k + 1) * Sz(i + 1, j + 1, k + 1)) /
                          (denom);

    // trilinear coefficients: f(x,y,z) = a+bx+cy+dz+exy+fxz+gyz+hxyz
    // Evaluate the coefficients
    a = CubeValues(0, 0, 0);
    b = CubeValues(1, 0, 0) - a;
    c = CubeValues(0, 1, 0) - a;
    d = CubeValues(0, 0, 1) - a;
    e = CubeValues(1, 1, 0) - a - b - c;
    f = CubeValues(1, 0, 1) - a - b - d;
    g = CubeValues(0, 1, 1) - a - c - d;
    h = CubeValues(1, 1, 1) - a - b - c - d - e - f - g;

    for (p = 0; p < npts; p++) {
        A = Points(p);
        x = A.x - 1.0 * i;
        y = A.y - 1.0 * j;
        z = A.z - 1.0 * k;
        CurveValues(p) = a + b * x + c * y + d * z + e * x * y + f * x * z +
                         g * y * z + h * x * y * z;
        //		printf("grad F = %f \n", sqrt(pow(Fx(i,j,k),2)+pow(Fy(i,j,k),2)+pow(Fz(i,j,k),2)));
    }

    integral = 0.0;

    for (p = 0; p < npts - 1; p++) {
        // Extract the line segment
        A = Points(p);
        B = Points(p + 1);
        vA = CurveValues(p);
        vB = CurveValues(p + 1);
        // Compute the length of the segment
        length = sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y) +
                      (A.z - B.z) * (A.z - B.z));
        integral += 0.5 * length * (vA + vB);
    }

    return integral;
}
//--------------------------------------------------------------------------------------------------------
inline void
pmmc_CubeSurfaceInterpVector(DoubleArray &Vec_x, DoubleArray &Vec_y,
                             DoubleArray &Vec_z, DoubleArray &CubeValues,
                             DTMutableList<Point> &Points, IntArray &Triangles,
                             DoubleArray &SurfaceVector, DoubleArray &VecAvg,
                             int i, int j, int k, int npts, int ntris) {
    Point A, B, C;
    int p;
    double vA, vB, vC;
    double x, y, z;
    double s, s1, s2, s3, temp;
    double a, b, c, d, e, f, g, h;

    // ................x component .............................
    // Copy the curvature values for the cube
    CubeValues(0, 0, 0) = Vec_x(i, j, k);
    CubeValues(1, 0, 0) = Vec_x(i + 1, j, k);
    CubeValues(0, 1, 0) = Vec_x(i, j + 1, k);
    CubeValues(1, 1, 0) = Vec_x(i + 1, j + 1, k);
    CubeValues(0, 0, 1) = Vec_x(i, j, k + 1);
    CubeValues(1, 0, 1) = Vec_x(i + 1, j, k + 1);
    CubeValues(0, 1, 1) = Vec_x(i, j + 1, k + 1);
    CubeValues(1, 1, 1) = Vec_x(i + 1, j + 1, k + 1);

    // trilinear coefficients: f(x,y,z) = a+bx+cy+dz+exy+fxz+gyz+hxyz
    a = CubeValues(0, 0, 0);
    b = CubeValues(1, 0, 0) - a;
    c = CubeValues(0, 1, 0) - a;
    d = CubeValues(0, 0, 1) - a;
    e = CubeValues(1, 1, 0) - a - b - c;
    f = CubeValues(1, 0, 1) - a - b - d;
    g = CubeValues(0, 1, 1) - a - c - d;
    h = CubeValues(1, 1, 1) - a - b - c - d - e - f - g;

    for (p = 0; p < npts; p++) {
        A = Points(p);
        x = A.x - 1.0 * i;
        y = A.y - 1.0 * j;
        z = A.z - 1.0 * k;
        SurfaceVector(p) = a + b * x + c * y + d * z + e * x * y + f * x * z +
                           g * y * z + h * x * y * z;
    }

    // ................y component .............................
    // Copy the curvature values for the cube
    CubeValues(0, 0, 0) = Vec_y(i, j, k);
    CubeValues(1, 0, 0) = Vec_y(i + 1, j, k);
    CubeValues(0, 1, 0) = Vec_y(i, j + 1, k);
    CubeValues(1, 1, 0) = Vec_y(i + 1, j + 1, k);
    CubeValues(0, 0, 1) = Vec_y(i, j, k + 1);
    CubeValues(1, 0, 1) = Vec_y(i + 1, j, k + 1);
    CubeValues(0, 1, 1) = Vec_y(i, j + 1, k + 1);
    CubeValues(1, 1, 1) = Vec_y(i + 1, j + 1, k + 1);

    // trilinear coefficients: f(x,y,z) = a+bx+cy+dz+exy+fxz+gyz+hxyz
    a = CubeValues(0, 0, 0);
    b = CubeValues(1, 0, 0) - a;
    c = CubeValues(0, 1, 0) - a;
    d = CubeValues(0, 0, 1) - a;
    e = CubeValues(1, 1, 0) - a - b - c;
    f = CubeValues(1, 0, 1) - a - b - d;
    g = CubeValues(0, 1, 1) - a - c - d;
    h = CubeValues(1, 1, 1) - a - b - c - d - e - f - g;

    for (p = 0; p < npts; p++) {
        A = Points(p);
        x = A.x - 1.0 * i;
        y = A.y - 1.0 * j;
        z = A.z - 1.0 * k;
        SurfaceVector(npts + p) = a + b * x + c * y + d * z + e * x * y +
                                  f * x * z + g * y * z + h * x * y * z;
    }

    // ................z component .............................
    // Copy the curvature values for the cube
    CubeValues(0, 0, 0) = Vec_z(i, j, k);
    CubeValues(1, 0, 0) = Vec_z(i + 1, j, k);
    CubeValues(0, 1, 0) = Vec_z(i, j + 1, k);
    CubeValues(1, 1, 0) = Vec_z(i + 1, j + 1, k);
    CubeValues(0, 0, 1) = Vec_z(i, j, k + 1);
    CubeValues(1, 0, 1) = Vec_z(i + 1, j, k + 1);
    CubeValues(0, 1, 1) = Vec_z(i, j + 1, k + 1);
    CubeValues(1, 1, 1) = Vec_z(i + 1, j + 1, k + 1);

    // trilinear coefficients: f(x,y,z) = a+bx+cy+dz+exy+fxz+gyz+hxyz
    a = CubeValues(0, 0, 0);
    b = CubeValues(1, 0, 0) - a;
    c = CubeValues(0, 1, 0) - a;
    d = CubeValues(0, 0, 1) - a;
    e = CubeValues(1, 1, 0) - a - b - c;
    f = CubeValues(1, 0, 1) - a - b - d;
    g = CubeValues(0, 1, 1) - a - c - d;
    h = CubeValues(1, 1, 1) - a - b - c - d - e - f - g;

    for (p = 0; p < npts; p++) {
        A = Points(p);
        x = A.x - 1.0 * i;
        y = A.y - 1.0 * j;
        z = A.z - 1.0 * k;
        SurfaceVector(2 * npts + p) = a + b * x + c * y + d * z + e * x * y +
                                      f * x * z + g * y * z + h * x * y * z;
    }
    //.............................................................................

    for (int r = 0; r < ntris; r++) {
        A = Points(Triangles(0, r));
        B = Points(Triangles(1, r));
        C = Points(Triangles(2, r));
        // Compute length of sides (assume dx=dy=dz)
        s1 = sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y) +
                  (A.z - B.z) * (A.z - B.z));
        s2 = sqrt((A.x - C.x) * (A.x - C.x) + (A.y - C.y) * (A.y - C.y) +
                  (A.z - C.z) * (A.z - C.z));
        s3 = sqrt((B.x - C.x) * (B.x - C.x) + (B.y - C.y) * (B.y - C.y) +
                  (B.z - C.z) * (B.z - C.z));
        s = 0.5 * (s1 + s2 + s3);
        temp = s * (s - s1) * (s - s2) * (s - s3);
        if (temp > 0.0) {
            // Increment the averaged values
            // x component
            vA = SurfaceVector(Triangles(0, r));
            vB = SurfaceVector(Triangles(1, r));
            vC = SurfaceVector(Triangles(2, r));
            VecAvg(0) += sqrt(temp) * 0.33333333333333333 * (vA + vB + vC);
            // y component
            vA = SurfaceVector(npts + Triangles(0, r));
            vB = SurfaceVector(npts + Triangles(1, r));
            vC = SurfaceVector(npts + Triangles(2, r));
            VecAvg(1) += sqrt(temp) * 0.33333333333333333 * (vA + vB + vC);
            // z component
            vA = SurfaceVector(2 * npts + Triangles(0, r));
            vB = SurfaceVector(2 * npts + Triangles(1, r));
            vC = SurfaceVector(2 * npts + Triangles(2, r));
            VecAvg(2) += sqrt(temp) * 0.33333333333333333 * (vA + vB + vC);
        }
    }
}
//--------------------------------------------------------------------------------------------------------
inline double pmmc_CubeSurfaceOrientation(DoubleArray &Orientation,
                                          DTMutableList<Point> &Points,
                                          IntArray &Triangles, int ntris) {
    int r;
    double temp, area, s, s1, s2, s3;
    double nx, ny, nz, normsq;
    Point A, B, C;
    area = 0.0;
    for (r = 0; r < ntris; r++) {
        A = Points(Triangles(0, r));
        B = Points(Triangles(1, r));
        C = Points(Triangles(2, r));
        // Compute the triangle normal vector
        nx = (B.y - A.y) * (C.z - A.z) - (B.z - A.z) * (C.y - A.y);
        ny = (B.z - A.z) * (C.x - A.x) - (B.x - A.x) * (C.z - A.z);
        nz = (B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x);
        normsq = 1.0 / (nx * nx + ny * ny + nz * nz);
        // Compute length of sides (assume dx=dy=dz)
        s1 = sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y) +
                  (A.z - B.z) * (A.z - B.z));
        s2 = sqrt((A.x - C.x) * (A.x - C.x) + (A.y - C.y) * (A.y - C.y) +
                  (A.z - C.z) * (A.z - C.z));
        s3 = sqrt((B.x - C.x) * (B.x - C.x) + (B.y - C.y) * (B.y - C.y) +
                  (B.z - C.z) * (B.z - C.z));
        s = 0.5 * (s1 + s2 + s3);
        temp = s * (s - s1) * (s - s2) * (s - s3);
        if (temp > 0.0) {
            temp = sqrt(temp);
            area += temp;
            Orientation(0) += temp * nx * nx * normsq; // Gxx
            Orientation(1) += temp * ny * ny * normsq; // Gyy
            Orientation(2) += temp * nz * nz * normsq; // Gzz
            Orientation(3) += temp * nx * ny * normsq; // Gxy
            Orientation(4) += temp * nx * nz * normsq; // Gxz
            Orientation(5) += temp * ny * nz * normsq; // Gyz
        }
    }
    return area;
}
//--------------------------------------------------------------------------------------------------------
inline double pmmc_CommonCurveSpeed(DoubleArray &CubeValues, DoubleArray &dPdt,
                                    DoubleArray &ReturnVector, DoubleArray &Fx,
                                    DoubleArray &Fy, DoubleArray &Fz,
                                    DoubleArray &Sx, DoubleArray &Sy,
                                    DoubleArray &Sz,
                                    DTMutableList<Point> &Points, int i, int j,
                                    int k, int npts) {
    NULL_USE(CubeValues);

    int p;
    double s, lwns, norm;
    double zeta;
    double tangent_x, tangent_y, tangent_z;
    double ns_x, ns_y, ns_z;
    double nwn_x, nwn_y, nwn_z;
    double nwns_x, nwns_y, nwns_z;
    Point P, A, B;
    lwns = 0.0;
    double ReturnValue = 0;
    NULL_USE(lwns);

    TriLinPoly Px, Py, Pz, SDx, SDy, SDz, Pt;
    Px.assign(Fx, i, j, k);
    Py.assign(Fy, i, j, k);
    Pz.assign(Fz, i, j, k);
    SDx.assign(Sx, i, j, k);
    SDy.assign(Sy, i, j, k);
    SDz.assign(Sz, i, j, k);
    Pt.assign(dPdt, i, j, k);

    for (p = 0; p < npts - 1; p++) {
        // Extract the line segment
        A = Points(p);
        B = Points(p + 1);
        // Midpoint of the line
        P.x = 0.5 * (A.x + B.x);
        P.y = 0.5 * (A.y + B.y);
        P.z = 0.5 * (A.z + B.z);
        // Compute the curve tangent
        tangent_x = A.x - B.x;
        tangent_y = A.y - B.y;
        tangent_z = A.z - B.z;

        // Get the normal to the solid surface
        ns_x = SDx.eval(P);
        ns_y = SDy.eval(P);
        ns_z = SDz.eval(P);
        norm = ns_x * ns_x + ns_y * ns_y + ns_z * ns_z;
        if (norm > 0.0) {
            ns_x /= norm;
            ns_y /= norm;
            ns_z /= norm;
        }

        // Get the Color gradient
        nwn_x = Px.eval(P);
        nwn_y = Py.eval(P);
        nwn_z = Pz.eval(P);
        // to compute zeta, consider the change only in the plane of the solid surface
        norm = nwn_x * ns_x + nwn_y * ns_y +
               nwn_z * ns_z; // component of color gradient aligned with solid
        nwns_x = nwn_x - norm * ns_x;
        nwns_y = nwn_y - norm * ns_y;
        nwns_z = nwn_z - norm * ns_z;
        // now {nwns_x, nwns_y, nwns_z} is the color gradient confined to the solid plane
        norm = sqrt(nwns_x * nwns_x + nwns_y * nwns_y + nwns_z * nwns_z);
        zeta = -Pt.eval(P) / norm;
        // normalize the normal to the common curve within the solid surface
        if (norm > 0.0) {
            nwns_x /= norm;
            nwns_y /= norm;
            nwns_z /= norm;
        }

        /*
		// normal to the wn interface
		norm = nwn_x*nwn_x + nwn_y*nwn_y + nwn_z*nwn_z;
		// Compute the interface speed
		if (norm > 0.0){
			nwn_x /= norm;
			nwn_y /= norm;
			nwn_z /= norm;
		}

		// Compute the normal to the common curve (vector product of tangent and solid normal)
		nwns_x = ns_y*tangent_z - ns_z*tangent_y;
		nwns_y = ns_z*tangent_x - ns_x*tangent_z;
		nwns_z = ns_x*tangent_y - ns_y*tangent_x;
		// dot product of nwns and nwn should be positive
		if (nwn_x*nwns_x+nwn_y*nwns_y+nwn_z*nwns_z < 0.0){
			nwns_x = -nwns_x;
			nwns_y = -nwns_y;
			nwns_z = -nwns_z;
		}
		// normalize the common curve normal
		norm = sqrt(nwns_x*nwns_x + nwns_y*nwns_y + nwns_z*nwns_z);
		if (norm > 0.0){
			nwns_x /= norm;
			nwns_y /= norm;
			nwns_z /= norm;
		}
		*/
        // Compute the length of the segment
        s = sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y) +
                 (A.z - B.z) * (A.z - B.z));
        // Add the length to the common line
        //	lwns += s;
        // Compute the common curve velocity
        if (norm > 0.0) {
            ReturnVector(0) += zeta * nwns_x * s;
            ReturnVector(1) += zeta * nwns_y * s;
            ReturnVector(2) += zeta * nwns_z * s;
            ReturnValue += zeta * s;
        }
    }
    NULL_USE(tangent_x);
    NULL_USE(tangent_y);
    NULL_USE(tangent_z);

    return ReturnValue;
}
inline void pmmc_CurveOrientation(DoubleArray &Orientation,
                                  DTMutableList<Point> &Points, int npts, int i,
                                  int j, int k) {

    double twnsx, twnsy, twnsz, length; // tangent, norm

    Point P, A, B;

    for (int p = 0; p < npts - 1; p++) {
        // Extract the line segment
        A = Points(p);
        B = Points(p + 1);
        P.x = 0.5 * (A.x + B.x) - 1.0 * i;
        P.y = 0.5 * (A.y + B.y) - 1.0 * j;
        P.z = 0.5 * (A.z + B.z) - 1.0 * k;

        A.x -= 1.0 * i;
        A.y -= 1.0 * j;
        A.z -= 1.0 * k;
        B.x -= 1.0 * i;
        B.y -= 1.0 * j;
        B.z -= 1.0 * k;

        length = sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y) +
                      (A.z - B.z) * (A.z - B.z)); // tangent vector
        twnsx = (B.x - A.x) / length;
        twnsy = (B.y - A.y) / length;
        twnsz = (B.z - A.z) / length;

        Orientation(0) += (1.0 - twnsx * twnsx) * length; // Gxx
        Orientation(1) += (1.0 - twnsy * twnsy) * length; // Gyy
        Orientation(2) += (1.0 - twnsz * twnsz) * length; // Gzz
        Orientation(3) += (0.0 - twnsx * twnsy) * length; // Gxy
        Orientation(4) += (0.0 - twnsx * twnsz) * length; // Gxz
        Orientation(5) += (0.0 - twnsy * twnsz) * length; // Gyz
    }
}

inline void pmmc_CurveCurvature(DoubleArray &f, DoubleArray &s,
                                DoubleArray &f_x, DoubleArray &f_y,
                                DoubleArray &f_z, DoubleArray &s_x,
                                DoubleArray &s_y, DoubleArray &s_z,
                                DoubleArray &KN, DoubleArray &KG, double &KNavg,
                                double &KGavg, DTMutableList<Point> &Points,
                                int npts, int ic, int jc, int kc) {
    NULL_USE(f);
    NULL_USE(s);
    NULL_USE(KN);
    NULL_USE(KG);

    int p, i, j, k;
    double length;
    double fxx, fyy, fzz, fxy, fxz, fyz, fx, fy, fz;
    double sxx, syy, szz, sxy, sxz, syz, sx, sy, sz;
    //	double Axx,Axy,Axz,Ayx,Ayy,Ayz,Azx,Azy,Azz;
    //	double Tx[8],Ty[8],Tz[8];	// Tangent vector
    //	double Nx[8],Ny[8],Nz[8];	// Principle normal
    double twnsx, twnsy, twnsz, nwnsx, nwnsy, nwnsz,
        K; // tangent,normal and curvature
    double nsx, nsy, nsz, norm;
    double nwx, nwy, nwz;
    double nwsx, nwsy, nwsz;
    double nwnx, nwny, nwnz;

    Point P, A, B;
    // Local trilinear approximation for tangent and normal vector
    TriLinPoly Tx, Ty, Tz, Nx, Ny, Nz, Sx, Sy, Sz;

    // Loop over the cube and compute the derivatives
    for (k = kc; k < kc + 2; k++) {
        for (j = jc; j < jc + 2; j++) {
            for (i = ic; i < ic + 2; i++) {
                sx = s_x(i, j, k);
                sy = s_y(i, j, k);
                sz = s_z(i, j, k);
                fx = f_x(i, j, k);
                fy = f_y(i, j, k);
                fz = f_z(i, j, k);

                // Normal to fluid surface
                Nx.Corners(i - ic, j - jc, k - kc) = fx;
                Ny.Corners(i - ic, j - jc, k - kc) = fy;
                Nz.Corners(i - ic, j - jc, k - kc) = fz;

                // Normal to solid surface
                Sx.Corners(i - ic, j - jc, k - kc) = sx;
                Sy.Corners(i - ic, j - jc, k - kc) = sy;
                Sz.Corners(i - ic, j - jc, k - kc) = sz;

                // Compute the tangent vector
                norm = sqrt((sy * fz - sz * fy) * (sy * fz - sz * fy) +
                            (sz * fx - sx * fz) * (sz * fx - sx * fz) +
                            (sx * fy - sy * fx) * (sx * fy - sy * fx));
                if (norm == 0.0)
                    norm = 1.0;
                Tx.Corners(i - ic, j - jc, k - kc) = (sy * fz - sz * fy) / norm;
                Ty.Corners(i - ic, j - jc, k - kc) = (sz * fx - sx * fz) / norm;
                Tz.Corners(i - ic, j - jc, k - kc) = (sx * fy - sy * fx) / norm;
            }
        }
    }

    // Assign the tri-linear polynomial coefficients
    Sx.assign();
    Sy.assign();
    Sz.assign();
    Tx.assign();
    Ty.assign();
    Tz.assign();
    Nx.assign();
    Ny.assign();
    Nz.assign();

    double tax, tay, taz, tbx, tby, tbz;

    for (p = 0; p < npts - 1; p++) {
        // Extract the line segment
        A = Points(p);
        B = Points(p + 1);
        P.x = 0.5 * (A.x + B.x) - 1.0 * i;
        P.y = 0.5 * (A.y + B.y) - 1.0 * j;
        P.z = 0.5 * (A.z + B.z) - 1.0 * k;

        A.x -= 1.0 * i;
        A.y -= 1.0 * j;
        A.z -= 1.0 * k;
        B.x -= 1.0 * i;
        B.y -= 1.0 * j;
        B.z -= 1.0 * k;

        length = sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y) +
                      (A.z - B.z) * (A.z - B.z));

        // tangent vector
        twnsx = B.x - A.x;
        twnsy = B.y - A.y;
        twnsz = B.z - A.z;
        norm = sqrt(twnsx * twnsx + twnsy * twnsy + twnsz * twnsz);
        if (norm == 0.0)
            norm = 1.0;
        twnsx /= norm;
        twnsy /= norm;
        twnsz /= norm;

        // *****************************************************
        // Compute tangent vector at A and B and normalize
        tax = Tx.eval(A);
        tay = Ty.eval(A);
        taz = Tz.eval(A);
        norm = sqrt(tax * tax + tay * tay + taz * taz);
        if (norm == 0.0)
            norm = 1.0;
        tax /= norm;
        tay /= norm;
        taz /= norm;

        tbx = Tx.eval(B);
        tby = Ty.eval(B);
        tbz = Tz.eval(B);
        norm = sqrt(tbx * tbx + tby * tby + tbz * tbz);
        if (norm == 0.0)
            norm = 1.0;
        tbx /= norm;
        tby /= norm;
        tbz /= norm;
        // *****************************************************

        // *****************************************************
        // Compute the divergence of the tangent vector
        nwnsx = (tax - tbx) / length;
        nwnsy = (tay - tby) / length;
        nwnsz = (taz - tbz) / length;
        K = sqrt(nwnsx * nwnsx + nwnsy * nwnsy + nwnsz * nwnsz);
        if (K == 0.0)
            K = 1.0;
        nwnsx /= K;
        nwnsy /= K;
        nwnsz /= K;
        // *****************************************************

        // Normal vector to the solid surface
        nsx = Sx.eval(P);
        nsy = Sy.eval(P);
        nsz = Sz.eval(P);
        norm = sqrt(nsx * nsx + nsy * nsy + nsz * nsz);
        if (norm == 0.0)
            norm = 1.0;
        nsx /= norm;
        nsy /= norm;
        nsz /= norm;

        // Normal vector to the fluid surface
        nwx = Nx.eval(P);
        nwy = Ny.eval(P);
        nwz = Nz.eval(P);
        norm = sqrt(nwx * nwx + nwy * nwy + nwz * nwz);
        if (norm == 0.0)
            norm = 1.0;
        nwx /= norm;
        nwy /= norm;
        nwz /= norm;

        // common curve normal in the solid surface tangent plane (rel. geodesic curvature)
        nwsx = twnsy * nsz - twnsz * nsy;
        nwsy = twnsz * nsx - twnsx * nsz;
        nwsz = twnsx * nsy - twnsy * nsx;
        norm = sqrt(nwsx * nwsx + nwsy * nwsy + nwsz * nwsz);
        if (norm == 0.0)
            norm = 1.0;
        nwsx /= norm;
        nwsy /= norm;
        nwsz /= norm;
        /* normal to ws interface boundary should point into fluid (same direction as gradient) */
        if (nwx * nwsx + nwy * nwsy + nwz * nwsz < 0.0) {
            nwsx = -nwsx;
            nwsy = -nwsy;
            nwsz = -nwsz;
        }

        // common curve normal in the fluid surface tangent plane (rel. geodesic curvature)
        nwnx = twnsy * nwz - twnsz * nwy;
        nwny = twnsz * nwx - twnsx * nwz;
        nwnz = twnsx * nwy - twnsy * nwx;
        norm = sqrt(nwnx * nwnx + nwny * nwny + nwnz * nwnz);
        if (norm == 0.0)
            norm = 1.0;
        nwnx /= norm;
        nwny /= norm;
        nwnz /= norm;
        /* normal to wn interface boundary should point into the solid */
        if (nsx * nwnx + nsy * nwny + nsz * nwnz > 0.0) {
            nwnx = -nwnx;
            nwny = -nwny;
            nwnz = -nwnz;
        }

        if (length > 0.0) {
            // normal curvature component in the direction of the solid surface
            //KNavg += K * (nsx * nwnsx + nsy * nwnsy + nsz * nwnsz) * length;
            KNavg += K * (nwnx * nwnsx + nwny * nwnsy + nwnz * nwnsz) * length;
            //geodesic curvature
            KGavg += K * (nwsx * nwnsx + nwsy * nwnsy + nwsz * nwnsz) * length;
        }
    }
    NULL_USE(fxx);
    NULL_USE(fyy);
    NULL_USE(fzz);
    NULL_USE(fxy);
    NULL_USE(fxz);
    NULL_USE(fyz);
    NULL_USE(fx);
    NULL_USE(fy);
    NULL_USE(fz);
    NULL_USE(sxx);
    NULL_USE(syy);
    NULL_USE(szz);
    NULL_USE(sxy);
    NULL_USE(sxz);
    NULL_USE(syz);
    NULL_USE(sx);
    NULL_USE(sy);
    NULL_USE(sz);
}

//--------------------------------------------------------------------------------------------------------
inline double pmmc_InterfaceSpeed(
    DoubleArray &dPdt, DoubleArray &P_x, DoubleArray &P_y, DoubleArray &P_z,
    DoubleArray &CubeValues, DTMutableList<Point> &Points, IntArray &Triangles,
    DoubleArray &SurfaceVector, DoubleArray &AvgSpeed, DoubleArray &AvgVel,
    int i, int j, int k, int npts, int ntris) {
    NULL_USE(CubeValues);
    NULL_USE(SurfaceVector);
    NULL_USE(npts);

    Point A, B, C, P;
    double x, y, z;
    double s, s1, s2, s3, area;
    double norm, zeta;
    double ReturnValue = 0.0;

    TriLinPoly Px, Py, Pz, Pt;
    Px.assign(P_x, i, j, k);
    Py.assign(P_y, i, j, k);
    Pz.assign(P_z, i, j, k);
    Pt.assign(dPdt, i, j, k);

    //.............................................................................
    // Compute the average speed of the interface
    for (int r = 0; r < ntris; r++) {
        A = Points(Triangles(0, r));
        B = Points(Triangles(1, r));
        C = Points(Triangles(2, r));
        // Compute length of sides (assume dx=dy=dz)
        s1 = sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y) +
                  (A.z - B.z) * (A.z - B.z));
        s2 = sqrt((A.x - C.x) * (A.x - C.x) + (A.y - C.y) * (A.y - C.y) +
                  (A.z - C.z) * (A.z - C.z));
        s3 = sqrt((B.x - C.x) * (B.x - C.x) + (B.y - C.y) * (B.y - C.y) +
                  (B.z - C.z) * (B.z - C.z));
        s = 0.5 * (s1 + s2 + s3);
        area = sqrt(s * (s - s1) * (s - s2) * (s - s3));
        // Compute the centroid P
        P.x = 0.33333333333333333 * (A.x + B.x + C.x);
        P.y = 0.33333333333333333 * (A.y + B.y + C.y);
        P.z = 0.33333333333333333 * (A.z + B.z + C.z);
        if (area > 0.0) {
            x = Px.eval(P);
            y = Py.eval(P);
            z = Pz.eval(P);
            norm = sqrt(x * x + y * y + z * z);
            if (norm == 0.0)
                norm = 1.0;
            // Compute the interface speed from time derivative and gradient (Level Set Equation)
            zeta = -Pt.eval(P) / norm;
            // Normalize the normal vector
            x /= norm;
            y /= norm;
            z /= norm;

            // Compute the average
            AvgVel(0) += area * zeta * x;
            AvgVel(1) += area * zeta * y;
            AvgVel(2) += area * zeta * z;
            AvgSpeed(r) = zeta * area;
            ReturnValue += zeta * area;
        }
    }
    return ReturnValue;
    //.............................................................................
}
//--------------------------------------------------------------------------------------------------------
inline double geomavg_EulerCharacteristic(DTMutableList<Point> &Points,
                                          IntArray &Triangles, int &npts,
                                          int &ntris, int &i, int &j, int &k) {
    NULL_USE(Points);
    NULL_USE(Triangles);
    NULL_USE(npts);
    NULL_USE(ntris);
    NULL_USE(i);
    NULL_USE(j);
    NULL_USE(k);

    /* REFERENCE
	 *  Huang, Liu, Lee, Yang, Tsang
	 *  On concise 3-D simple point comparisons: a marching cubes paradigm
	 *  IEEE Transactions on Medical Imaging, Vol. 28, No. 1
	 *  January 2009
	 */
    // Compute the Euler characteristic for triangles in a cube
    /*
	bool graph[3][3];

	double CountSideExternal=0;
	double ShareSideInternal=0;

	// Exclude edges and vertices shared with between multiple cubes
	for (int p=0; p<ntris; p++){
		Point A = Points(Triangles(0,p));
		Point B = Points(Triangles(1,p));
		Point C = Points(Triangles(2,p));

		// Check to see if triangle sides are on the cube edge
		if ( fabs(floor(A.x) - A.x) < 1e-14 && fabs(floor(B.x) - B.x) < 1e-14 )		CountSideExternal+=1.0;
		if ( fabs(floor(A.x) - A.x) < 1e-14 && fabs(floor(C.x) - C.x) < 1e-14 )		CountSideExternal+=1.0;
		if ( fabs(floor(B.x) - B.x) < 1e-14 && fabs(floor(C.x) - C.x) < 1e-14 )		CountSideExternal+=1.0;

		// See if any edges were already counted
		for (int r=0; r<p; r++){
			Point D = Points(Triangles(0,r));
			Point E = Points(Triangles(1,r));
			Point F = Points(Triangles(2,r));

			// Initialize the graph
			for (int jj=0; jj<3; jj++){
				for (int ii=0; ii<3; ii++){
					graph[ii][jj] = false;
				}
			}

			// Check to see if the triangles share vertices
			if (Triangles(0,r) == Triangles(0,p))	graph[0][0]=true;
			if (Triangles(0,r) == Triangles(1,p))	graph[0][1]=true;
			if (Triangles(0,r) == Triangles(2,p))	graph[0][2]=true;

			if (Triangles(1,r) == Triangles(0,p))	graph[1][0]=true;
			if (Triangles(1,r) == Triangles(1,p))	graph[1][1]=true;
			if (Triangles(1,r) == Triangles(2,p))	graph[1][2]=true;

			if (Triangles(2,r) == Triangles(0,p))	graph[2][0]=true;
			if (Triangles(2,r) == Triangles(1,p))	graph[2][1]=true;
			if (Triangles(2,r) == Triangles(2,p))	graph[2][2]=true;

			// Count the number of vertices that are shared
			int count=0;
			for (int jj=0; jj<3; jj++){
				for (int ii=0; ii<3; ii++){
					if (graph[ii][jj] == true) count++;
				}
			}
			// 0,1,2 are valid possible numbers for shared vertices
			// if 3 are obtained, then it is the same triangle :(
			if (count == 2)  ShareSideInternal += 1.0;
		}
	}
	*/
    double EulerChar, nvert, nface;

    // Number of faces is number of triangles
    nface = double(ntris);
    // Each vertex is shared by four cubes
    nvert = double(npts);
    // Subtract shared sides to avoid double counting
    //	nside = 2.0*double(npts)- 3.0 - 0.5*double(npts);
    double nside_extern = double(npts);
    double nside_intern = double(npts) - 3.0;
    EulerChar = 0.0;
    if (npts > 0)
        EulerChar = (0.25 * nvert - nside_intern - 0.5 * nside_extern + nface);
    return EulerChar;
}

inline double mink_phase_epc6(IntArray &PhaseID, DoubleArray &CubeValues,
                              int PhaseLabel, int &i, int &j, int &k) {

    /*
	 *  Compute the Euler Poincare characteristic for a phase based on 6 adjacency
	 *  compute the local contribution within the specified cube must sum to get total
	 *  	PhaseID -- label of phases on a mesh
	 *  	PhaseLabel -- label assigned to phase that we are computing
	 *  	i,j,k -- identifies the local cube
	 */

    double epc_ncubes, epc_nface, epc_nedge, epc_nvert;
    epc_ncubes = epc_nface = epc_nedge = epc_nvert = 0;

    // Assign CubeValues to be phase indicator for
    for (int kk = 0; kk < 2; kk++) {
        for (int jj = 0; jj < 2; jj++) {
            for (int ii = 0; ii < 2; ii++) {
                if (PhaseID(i + ii, j + jj, k + kk) == PhaseLabel)
                    CubeValues(ii, jj, kk) = 1;
                else
                    CubeValues(ii, jj, kk) = 0;
            }
        }
    }
    // update faces edges and cubes for NWP
    epc_ncubes += CubeValues(0, 0, 0) * CubeValues(1, 0, 0) *
                  CubeValues(0, 1, 0) * CubeValues(0, 0, 1) *
                  CubeValues(1, 1, 0) * CubeValues(1, 0, 1) *
                  CubeValues(0, 1, 1) * CubeValues(1, 1, 1);
    //  three faces (others shared by other cubes)
    epc_nface += CubeValues(1, 0, 0) * CubeValues(1, 1, 0) *
                 CubeValues(1, 0, 1) * CubeValues(1, 1, 1);
    epc_nface += CubeValues(0, 1, 0) * CubeValues(1, 1, 0) *
                 CubeValues(0, 1, 1) * CubeValues(1, 1, 1);
    epc_nface += CubeValues(0, 0, 1) * CubeValues(0, 1, 1) *
                 CubeValues(1, 0, 1) * CubeValues(1, 1, 1);
    // six of twelve edges (others shared by other cubes)
    epc_nedge += CubeValues(1, 1, 1) * CubeValues(1, 1, 0);
    epc_nedge += CubeValues(1, 1, 1) * CubeValues(1, 0, 1);
    epc_nedge += CubeValues(1, 1, 1) * CubeValues(0, 1, 1);
    epc_nedge += CubeValues(1, 0, 1) * CubeValues(1, 0, 0);
    epc_nedge += CubeValues(1, 0, 1) * CubeValues(0, 0, 1);
    epc_nedge += CubeValues(0, 1, 1) * CubeValues(0, 0, 1);
    // four of eight vertices
    epc_nvert += CubeValues(1, 1, 0);
    epc_nvert += CubeValues(1, 0, 1);
    epc_nvert += CubeValues(0, 1, 1);
    epc_nvert += CubeValues(1, 1, 1);

    double chi = epc_nvert - epc_nedge + epc_nface - epc_ncubes;
    return chi;
}

inline double mink_EulerCharacteristic(DTMutableList<Point> &Points,
                                       IntArray &Triangles,
                                       DoubleArray &CubeValues, int &npts,
                                       int &ntris, int &i, int &j, int &k) {

    NULL_USE(CubeValues);

    // Compute the Euler characteristic for triangles in a cube
    // Exclude edges and vertices shared with between multiple cubes
    double EulerChar;
    int nvert = npts;
    int nside = 2 * nvert - 3;
    int nface = nvert - 2;

    //if (ntris != nface){
    //	nface = ntris;
    //	nside =
    //}
    //...........................................................
    // Check that this point is not on a previously computed face
    // Note direction that the marching cubes algorithm marches
    // In parallel, other sub-domains fill in the lower boundary
    for (int p = 0; p < npts; p++) {
        Point PT = Points(p);
        if (PT.x - double(i) < 1e-12)
            nvert -= 1;
        else if (PT.y - double(j) < 1e-12)
            nvert -= 1;
        else if (PT.z - double(k) < 1e-12)
            nvert -= 1;
    }
    // Remove redundantly computed edges (shared by two cubes across a cube face)
    for (int p = 0; p < ntris; p++) {
        Point A = Points(Triangles(0, p));
        Point B = Points(Triangles(1, p));
        Point C = Points(Triangles(2, p));

        // Check side A-B
        bool newside = true;
        if (A.x - double(i) < 1e-12 && B.x - double(i) < 1e-12)
            newside = false;
        if (A.y - double(j) < 1e-12 && B.y - double(j) < 1e-12)
            newside = false;
        if (A.z - double(k) < 1e-12 && B.z - double(k) < 1e-12)
            newside = false;
        if (!newside)
            nside -= 1;

        // Check side A-C
        newside = true;
        if (A.x - double(i) < 1e-12 && C.x - double(i) < 1e-12)
            newside = false;
        if (A.y - double(j) < 1e-12 && C.y - double(j) < 1e-12)
            newside = false;
        if (A.z - double(k) < 1e-12 && C.z - double(k) < 1e-12)
            newside = false;
        if (!newside)
            nside -= 1;

        // Check side B-C
        newside = true;
        if (B.x - double(i) < 1e-12 && C.x - double(i) < 1e-12)
            newside = false;
        if (B.y - double(j) < 1e-12 && C.y - double(j) < 1e-12)
            newside = false;
        if (B.z - double(k) < 1e-12 && C.z - double(k) < 1e-12)
            newside = false;
        if (!newside)
            nside -= 1;
    }

    EulerChar = 1.0 * (nvert - nside + nface);
    return EulerChar;
}

#endif

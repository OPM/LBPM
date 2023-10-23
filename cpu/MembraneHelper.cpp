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

extern "C" void Membrane_D3Q19_Unpack(int q, int *list, int *links, int start,
                                      int linkCount, double *recvbuf,
                                      double *dist, int N) {
    //....................................................................................
    // Unack distribution from the recv buffer
    // Distribution q matche Cqx, Cqy, Cqz
    // swap rule means that the distributions in recvbuf are OPPOSITE of q
    // dist may be even or odd distributions stored by stream layout
    //....................................................................................
    int n, idx, link;
    for (link = 0; link < linkCount; link++) {

        idx = links[start + link];
        // Get the value from the list -- note that n is the index is from the send (non-local) process
        n = list[start + idx];
        // unpack the distribution to the proper location
        if (!(n < 0))
            dist[q * N + n] = recvbuf[start + idx];
    }
}

extern "C" void Membrane_D3Q19_Transport(int q, int *list, int *links,
                                         double *coef, int start, int offset,
                                         int linkCount, double *recvbuf,
                                         double *dist, int N) {
    //....................................................................................
    // Unack distribution from the recv buffer
    // Distribution q matche Cqx, Cqy, Cqz
    // swap rule means that the distributions in recvbuf are OPPOSITE of q
    // dist may be even or odd distributions stored by stream layout
    //....................................................................................
    int n, idx, link;
    double alpha;
    for (link = offset; link < linkCount; link++) {

        idx = list[start + link];
        // Get the value from the list -- note that n is the index is from the send (non-local) process
        n = list[start + idx];
        alpha = coef[start + idx];
        // unpack the distribution to the proper location
        if (!(n < 0))
            dist[q * N + n] = alpha * recvbuf[start + idx];
    }
}

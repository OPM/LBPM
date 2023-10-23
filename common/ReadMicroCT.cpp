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
#include "common/ReadMicroCT.h"
#include "common/Utilities.h"

#include "zlib.h"

#include <cstring>

// Read a file into memory
std::vector<char> readFile(const std::string &filename) {
    auto fid = fopen(filename.c_str(), "rb");
    INSIST(fid, "File does not exist: " + filename);
    fseek(fid, 0, SEEK_END);
    size_t bytes = ftell(fid);
    fseek(fid, 0, SEEK_SET);
    std::vector<char> data(bytes);
    size_t bytes2 = fread(data.data(), 1, bytes, fid);
    ASSERT(bytes == bytes2);
    fclose(fid);
    return data;
}

// Decompress a gzip buffer
std::vector<char> gunzip(const std::vector<char> &in) {
    z_stream stream;
    std::vector<char> out(1000000);
    stream.next_in = (Bytef *)in.data();
    stream.avail_in = in.size();
    stream.total_in = 0;
    stream.zalloc = Z_NULL;
    stream.zfree = Z_NULL;
    stream.opaque = Z_NULL;
    stream.next_out = (Bytef *)out.data();
    stream.avail_out = out.size();
    stream.total_out = 0;
    ASSERT(inflateInit2(&stream, 16 + MAX_WBITS) == Z_OK);
    bool finished = inflate(&stream, Z_SYNC_FLUSH) == Z_STREAM_END;
    while (!finished && stream.msg == Z_NULL) {
        out.resize(2 * out.size());
        stream.next_out = (Bytef *)&out[stream.total_out];
        stream.avail_out = out.size() - stream.total_out;
        finished = inflate(&stream, Z_SYNC_FLUSH) == Z_STREAM_END;
    }
    ASSERT(stream.msg == Z_NULL);
    out.resize(stream.total_out);
    inflateEnd(&stream);
    return out;
}

// Read the compressed micro CT data
Array<uint8_t> readMicroCT(const std::string &filename) {
    auto in = readFile(filename);
    auto out = gunzip(in);
    ASSERT(out.size() == 1024 * 1024 * 1024);
    Array<uint8_t> data(1024, 1024, 1024);
    memcpy(data.data(), out.data(), data.length());
    return data;
}

// Read the compressed micro CT data and distribute
Array<uint8_t> readMicroCT(const Database &domain, const Utilities::MPI &comm) {
    // Get the local problem info
    auto n = domain.getVector<int>("n");
    int rank = comm.getRank();
    auto nproc = domain.getVector<int>("nproc");
    RankInfoStruct rankInfo(rank, nproc[0], nproc[1], nproc[2]);

    // Determine the largest file number to get
    int Nfx = (n[0] * rankInfo.nx + 1023) / 1024;
    int Nfy = (n[1] * rankInfo.ny + 1023) / 1024;
    int Nfz = (n[2] * rankInfo.nz + 1023) / 1024;

    // Load one of the files if rank < largest file
    Array<uint8_t> data;
    RankInfoStruct srcRankInfo(rank, Nfx, Nfy, Nfz);
    if (srcRankInfo.ix >= 0) {
        auto filename = domain.getScalar<std::string>("Filename");
        char tmp[100];
        if (filename.find("0x_0y_0z.gbd.gz") != std::string::npos) {
            sprintf(tmp, "%ix_%iy_%iz.gbd.gz", srcRankInfo.ix, srcRankInfo.jy,
                    srcRankInfo.kz);
            filename = filename.replace(filename.find("0x_0y_0z.gbd.gz"), 15,
                                        std::string(tmp));
        } else if (filename.find("x0_y0_z0.gbd.gz") != std::string::npos) {
            sprintf(tmp, "x%i_y%i_z%i.gbd.gz", srcRankInfo.ix, srcRankInfo.jy,
                    srcRankInfo.kz);
            filename = filename.replace(filename.find("x0_y0_z0.gbd.gz"), 15,
                                        std::string(tmp));
        } else {
            ERROR("Invalid name for first file");
        }
        data = readMicroCT(filename);
    }

    // Redistribute the data
    data = redistribute(srcRankInfo, data, rankInfo, {n[0], n[1], n[2]}, comm);

    // Relabel the data
    auto ReadValues = domain.getVector<int>("ReadValues");
    auto WriteValues = domain.getVector<int>("WriteValues");
    ASSERT(ReadValues.size() == WriteValues.size());
    int readMaxValue = 0;
    for (auto v : ReadValues)
        readMaxValue = std::max(data.max() + 1, v);
    std::vector<int> map(readMaxValue + 1, -1);
    for (size_t i = 0; i < ReadValues.size(); i++)
        map[ReadValues[i]] = WriteValues[i];
    for (size_t i = 0; i < data.length(); i++) {
        int t = data(i);
        ASSERT(t >= 0 && t <= readMaxValue);
        data(i) = map[t];
    }

    return data;
}

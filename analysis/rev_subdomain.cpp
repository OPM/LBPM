#include "analysis/rev_analysis.h"
#include "analysis/rev_subdomain_flow.h"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>
#include <vector>

namespace {
using Shape = std::array<int, 3>;
constexpr double porosity_weight = 0.6;
constexpr double surface_weight = 0.4;
constexpr double permeability_weight = 0.6;
constexpr double tortuosity_weight = 0.4;

double relativeDeviation(double local, double target) {
    return target != 0.0 ? std::abs(local - target) / std::abs(target)
        : (local == 0.0 ? 0.0 : std::numeric_limits<double>::infinity());
}

struct Selection {
    Shape start{{0, 0, 0}};
    double value = 0.0;
    double deviation = std::numeric_limits<double>::infinity();
    double secondary_value = 0.0;
    double secondary_deviation = 0.0;
    double score = std::numeric_limits<double>::infinity();
    uint64_t candidates = 0;
    bool found = false;

    void considerHydraulics(const Shape &position, double local, double target,
                            double local_tortuosity, double target_tortuosity) {
        if (!std::isfinite(local))
            throw std::runtime_error("Subdomain selection received a nonfinite property");
        if (!std::isfinite(local_tortuosity)) {
            ++candidates; // Undefined tortuosity is evaluated but cannot be ranked.
            return;
        }
        const double relative = relativeDeviation(local, target);
        const double tortuosity_relative = relativeDeviation(local_tortuosity, target_tortuosity);
        update(position, local, relative, local_tortuosity, tortuosity_relative,
               permeability_weight * relative + tortuosity_weight * tortuosity_relative);
    }

    void considerGeometry(const Shape &position, double local, double target,
                          double local_surface, double target_surface) {
        if (!std::isfinite(local) || !std::isfinite(local_surface))
            throw std::runtime_error("Subdomain selection received a nonfinite geometric property");
        const double relative = relativeDeviation(local, target);
        const double surface_relative = relativeDeviation(local_surface, target_surface);
        update(position, local, relative, local_surface, surface_relative,
               porosity_weight * relative + surface_weight * surface_relative);
    }

    void update(const Shape &position, double local, double relative,
                double secondary, double secondary_relative, double candidate_score) {
        ++candidates;
        // Strict comparison preserves the first (z, y, x) position on ties.
        if (!found || candidate_score < score) {
            found = true;
            start = position;
            value = local;
            deviation = relative;
            secondary_value = secondary;
            secondary_deviation = secondary_relative;
            score = candidate_score;
        }
    }
};

std::string csvString(const std::string &value) {
    std::string result = "\"";
    for (char c : value) {
        result += c;
        if (c == '"') result += '"';
    }
    return result + '"';
}

std::string rawFilename(const char *property, const Shape &size) {
    return std::string("subdomain_") + property + '_' + std::to_string(size[0]) + 'x' +
        std::to_string(size[1]) + 'x' + std::to_string(size[2]) + ".raw";
}

void protectSource(const std::string &filename,
                   const std::array<std::string, 2> &raw_files) {
    // Check identity, including symlink and hard-link aliases, before any
    // result or staging file can be opened for writing.
    struct stat original;
    if (::stat(filename.c_str(), &original) != 0)
        return; // Missing sources are diagnosed when an export is needed.
    const auto check_output = [&](const std::string &reserved) {
        struct stat output;
        if (::stat(reserved.c_str(), &output) == 0 && original.st_dev == output.st_dev &&
            original.st_ino == output.st_ino)
            throw std::runtime_error("Selection output aliases the original image: " + reserved);
    };
    for (const char *reserved : {"subdomain_selection.csv", "geometrical-rev-domain.db",
                                 "hydraulic-rev-domain.db", "geometrical-rev-domain.db.tmp",
                                 "hydraulic-rev-domain.db.tmp"}) {
        check_output(reserved);
    }
    for (const auto &raw_file : raw_files) {
        if (raw_file.empty()) continue;
        check_output(raw_file);
        check_output(raw_file + ".tmp");
    }
}

// Domain::Decomp maps labels (and narrows 16-bit images). Export from the
// original file instead so selection never changes the source representation.
struct RawSource {
    std::string filename, read_type;
    Shape shape, offset{{0, 0, 0}};
    int bytes = 1;
    std::ifstream input;

    RawSource(const Database &domain, const Shape &interior) {
        if (domain.keyExists("Filename")) {
            filename = domain.getScalar<std::string>("Filename");
            read_type = domain.getWithDefault<std::string>("ReadType", "8bit");
            if (read_type != "8bit" && read_type != "16bit")
                throw std::runtime_error("Subdomain export supports raw 8bit/16bit geometry");
            bytes = read_type == "16bit" ? 2 : 1;
            auto n = domain.getVector<int>("N");
            if (n.size() != 3)
                throw std::runtime_error("Subdomain export requires three Domain.N dimensions");
            std::copy(n.begin(), n.end(), shape.begin());
            if (domain.keyExists("offset")) {
                auto origin = domain.getVector<int>("offset");
                if (origin.size() != 3)
                    throw std::runtime_error("Domain.offset must have three coordinates");
                std::copy(origin.begin(), origin.end(), offset.begin());
            }
        } else if (domain.keyExists("GridFile")) {
            throw std::runtime_error("Subdomain export requires raw geometry or ID.00000; GridFile is unsupported");
        } else {
            if (domain.getWithDefault<int>("BC", 0) != 0)
                throw std::runtime_error("Subdomain selection of ID.00000 requires BC = 0 to avoid implicit inlet/outlet layers");
            filename = "ID.00000";
            read_type = "signed8bit";
            for (int d = 0; d < 3; ++d) {
                shape[d] = interior[d] + 2;
                offset[d] = 1;
            }
        }
        for (int d = 0; d < 3; ++d)
            if (shape[d] <= 0 || offset[d] < 0 ||
                int64_t(offset[d]) + interior[d] > shape[d])
                throw std::runtime_error("Analyzed domain exceeds the original raw image bounds");
        input.open(filename, std::ios::binary | std::ios::ate);
        const auto expected = std::streamoff(shape[0]) * shape[1] * shape[2] * bytes;
        if (!input || input.tellg() != expected)
            throw std::runtime_error("Raw source size does not match its dimensions and data type: " + filename);
    }

    Array<signed char> poreMask(const Database &domain, const Shape &interior) {
        // Decomp changes inlet/outlet voxels after relabeling. Recover the
        // pre-modification pore mask using exactly its signed-byte convention,
        // including the narrowing of 16-bit labels before label mapping.
        std::array<signed char, 256> labels;
        for (size_t i = 0; i < labels.size(); ++i)
            labels[i] = static_cast<signed char>(i);
        if (domain.keyExists("Filename")) {
            const auto read = domain.getVector<int>("ReadValues");
            const auto write = domain.getVector<int>("WriteValues");
            if (read.size() != write.size())
                throw std::runtime_error("ReadValues and WriteValues must have equal lengths");
            std::array<bool, 256> mapped{};
            for (size_t i = 0; i < read.size(); ++i) {
                const auto index = static_cast<unsigned char>(read[i]);
                if (!mapped[index]) {
                    labels[index] = static_cast<signed char>(write[i]);
                    mapped[index] = true;
                }
            }
        }
        Array<signed char> pores(interior[0], interior[1], interior[2]);
        std::vector<char> row(size_t(interior[0]) * bytes);
        for (int z = 0; z < interior[2]; ++z) {
            for (int y = 0; y < interior[1]; ++y) {
                const auto index = (std::streamoff(z + offset[2]) * shape[1] +
                                    y + offset[1]) * shape[0] + offset[0];
                input.seekg(index * bytes);
                input.read(row.data(), std::streamsize(row.size()));
                if (!input)
                    throw std::runtime_error("Failed to read original rock pore labels");
                for (int x = 0; x < interior[0]; ++x) {
                    unsigned char label;
                    if (bytes == 2) {
                        uint16_t original;
                        std::memcpy(&original, row.data() + size_t(x) * 2, 2);
                        label = static_cast<unsigned char>(original);
                    } else {
                        label = static_cast<unsigned char>(row[x]);
                    }
                    pores(x, y, z) = labels[label] > 0 ? 1 : 0;
                }
            }
        }
        return pores;
    }

    void save(const std::string &output, const Shape &start, const Shape &size) {
        const std::string temporary = output + ".tmp";
        std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
        if (!out) throw std::runtime_error("Cannot write selected raw image: " + temporary);
        std::vector<char> row(size_t(size[0]) * bytes);
        for (int z = 0; z < size[2]; ++z) {
            for (int y = 0; y < size[1]; ++y) {
                auto index = (std::streamoff(start[2] + offset[2] + z) * shape[1] +
                              start[1] + offset[1] + y) * shape[0] + start[0] + offset[0];
                input.seekg(index * bytes);
                input.read(row.data(), std::streamsize(row.size()));
                out.write(row.data(), std::streamsize(row.size()));
                if (!input || !out)
                    throw std::runtime_error("Failed to copy original subdomain bytes");
            }
        }
        out.close();
        if (!out || std::rename(temporary.c_str(), output.c_str()) != 0)
            throw std::runtime_error("Cannot finish selected raw image: " + output);
    }

    void saveDomain(const std::string &output, const Database &original_domain,
                    const Shape &start, const Shape &size, double voxel_length) const {
        auto selected = original_domain.cloneDatabase();
        char *resolved = ::realpath(filename.c_str(), nullptr);
        if (!resolved)
            throw std::runtime_error("Cannot resolve original geometry path: " + filename);
        const std::string absolute_filename(resolved);
        std::free(resolved);
        selected->putScalar<std::string>("Filename", absolute_filename);
        selected->putScalar<std::string>("ReadType", bytes == 2 ? "16bit" : "8bit");
        selected->putVector<int>("N", {shape[0], shape[1], shape[2]});
        selected->putVector<int>("n", {size[0], size[1], size[2]});
        selected->putVector<int>("offset", {offset[0] + start[0], offset[1] + start[1],
                                             offset[2] + start[2]});
        selected->putVector<int>("nproc", {1, 1, 1});
        selected->putScalar<double>("voxel_length", voxel_length);
        if (selected->keyExists("L"))
            selected->putVector<double>("L", {size[0] * voxel_length, size[1] * voxel_length,
                                               size[2] * voxel_length});
        selected->putScalar<int>("BC", 0);
        selected->putVector<int>("InletLayers", {0, 0, 0});
        selected->putVector<int>("OutletLayers", {0, 0, 0});
        if (!original_domain.keyExists("Filename")) {
            // ID.00000 already stores mapped signed-byte labels. Reading it as
            // a raw full image must preserve that mapping and include its halo
            // in N and offset; unmatched labels also pass through unchanged.
            selected->putVector<int>("ReadValues", {0, 1, 2});
            selected->putVector<int>("WriteValues", {0, 1, 2});
        }
        Database result;
        result.putDatabase("Domain", selected);
        const std::string temporary = output + ".tmp";
        std::ofstream out(temporary, std::ios::trunc);
        if (!out)
            throw std::runtime_error("Cannot write selected Domain database: " + temporary);
        out << std::setprecision(std::numeric_limits<double>::max_digits10);
        result.print(out);
        out.close();
        if (!out || std::rename(temporary.c_str(), output.c_str()) != 0)
            throw std::runtime_error("Cannot finish selected Domain database: " + output);
    }
};
} // namespace

void ScaLBL_REV::SubDomainSelectImpl(const ModelData &model) {
    if (comm.getSize() != 1)
        throw std::runtime_error("Subdomain selection requires one MPI process");
    const Shape image{{model.Nx - 2, model.Ny - 2, model.Nz - 2}};
    const double h = model.voxel_length;
    if (!std::isfinite(h) || h <= 0.0)
        throw std::runtime_error("Subdomain selection requires a positive finite voxel length");
    auto domain = db->getDatabase("Domain");
    const auto configured_size = domain->getVector<int>("n");
    if (configured_size.size() != 3 ||
        !std::equal(configured_size.begin(), configured_size.end(), image.begin()))
        throw std::runtime_error("Selection model dimensions do not match Domain.n");
    if (subdomain_size[0] > 0)
        for (int d = 0; d < 3; ++d)
            if (subdomain_size[d] > image[d])
                throw std::runtime_error("REV.subdomain_size exceeds the analyzed Domain.n dimensions");
    const auto selection_dimensions = [&](const StatisticalREVResult &stat) {
        return subdomain_size[0] > 0 ? subdomain_size : stat.dimensions;
    };
    bool synthetic_layers = false;
    for (const char *key : {"InletLayers", "OutletLayers"}) {
        if (domain->keyExists(key)) {
            for (int layers : domain->getVector<int>(key))
                synthetic_layers = synthetic_layers || layers != 0;
        }
    }
    std::array<std::string, 2> raw_files;
    if (save_subdomain_raw) {
        if (geometrical && stat_poro_rev.size > 0)
            raw_files[0] = rawFilename("porosity", selection_dimensions(stat_poro_rev));
        if (hydraulic && stat_perm_rev.size > 0)
            raw_files[1] = rawFilename("permeability", selection_dimensions(stat_perm_rev));
    }
    std::unique_ptr<RawSource> source;
    protectSource(domain->getWithDefault<std::string>("Filename", "ID.00000"), raw_files);
    if ((geometrical && stat_poro_rev.size > 0) ||
        (hydraulic && stat_perm_rev.size > 0))
        source.reset(new RawSource(*domain, image));

    Array<signed char> original_pores;
    if (synthetic_layers && source) {
        original_pores = source->poreMask(*domain, image);
        std::cout << "Subdomain selection uses original rock labels before inlet/outlet modifications; "
                     "statistical REV targets are retained." << std::endl;
    }
    const auto is_pore = [&](int x, int y, int z) {
        return synthetic_layers ? original_pores(x, y, z) > 0
            : model.Distance(x + 1, y + 1, z + 1) > 0.0;
    };

    std::ofstream summary("subdomain_selection.csv", std::ios::trunc);
    if (!summary) throw std::runtime_error("Cannot write subdomain_selection.csv");
    summary << std::setprecision(std::numeric_limits<double>::max_digits10)
            << "property,status,rev_size_vx,rev_size_um,nx,ny,nz,lx_um,ly_um,lz_um,"
               "x0,y0,z0,x1,y1,z1,statistical_mean,selected_value,relative_deviation,"
               "candidate_count,raw_file,source_file,read_type,bytes_per_voxel,"
               "source_nx,source_ny,source_nz,source_x0,source_y0,source_z0,"
               "source_x1,source_y1,source_z1,domain_file,statistical_mean_surface,"
               "selected_surface,relative_deviation_surface,selection_score,porosity_weight,surface_weight,"
               "statistical_mean_tortuosity,selected_tortuosity,relative_deviation_tortuosity,"
               "permeability_weight,tortuosity_weight\n";

    const auto write_status = [&](const char *name, const char *status,
                                  const StatisticalREVResult *stat, uint64_t candidates) {
        summary << name << ',' << status;
        const Shape size = stat ? selection_dimensions(*stat) : Shape{{0, 0, 0}};
        for (int column = 2; column < 45; ++column) {
            summary << ',';
            if (column == 19) summary << candidates;
            if (!stat) continue;
            if (column == 2) summary << stat->size;
            else if (column == 3) summary << stat->size * h;
            else if (column >= 4 && column <= 6) summary << size[column - 4];
            else if (column >= 7 && column <= 9) summary << size[column - 7] * h;
            else if (column == 16) summary << stat->mean;
            else if (column == 40 && std::isfinite(stat->tortuosity_mean))
                summary << stat->tortuosity_mean;
        }
        summary << '\n';
        summary.flush();
        if (!summary) throw std::runtime_error("Failed to write subdomain_selection.csv");
    };

    bool device_ready = false;
    for (int property = 0; property < 2; ++property) {
        const bool porosity = property == 0;
        const char *name = porosity ? "porosity" : "permeability";
        const bool enabled = porosity ? geometrical : hydraulic;
        const auto &stat = porosity ? stat_poro_rev : stat_perm_rev;
        if (!enabled || stat.size == 0) {
            write_status(name, enabled ? "not_converged" : "disabled", nullptr, 0);
            continue;
        }
        const Shape size = selection_dimensions(stat);
        for (int d = 0; d < 3; ++d)
            if (size[d] < 1 || size[d] > image[d])
                throw std::runtime_error("Statistical REV dimensions exceed the input image");
        if (!std::isfinite(stat.mean))
            throw std::runtime_error("Statistical REV mean must be finite");
        if (!porosity && !std::isfinite(stat.tortuosity_mean)) {
            write_status(name, "unavailable_tortuosity", &stat, 0);
            std::cout << "Hydraulic selection has no valid statistical tortuosity target at the permeability REV size."
                      << std::endl;
            continue;
        }
        const double volume = double(size[0]) * size[1] * size[2];
        const uint64_t expected = uint64_t(image[0] - size[0] + 1) *
                                  (image[1] - size[1] + 1) * (image[2] - size[2] + 1);
        std::cout << "Selecting " << name << " subdomain: " << expected
                  << " windows of " << size[0] << 'x' << size[1] << 'x' << size[2] << std::endl;
        Selection best;
        if (porosity) {
            if (!std::isfinite(stat.surface_mean) || stat.surface_mean < 0.0)
                throw std::runtime_error("Geometric selection requires the statistical surface mean at the porosity REV size");
            // Four sliding z slabs and summed-area tables count pore voxels
            // and internal x/y/z interface faces, with O(Nx*Ny) memory.
            const size_t pitch = size_t(image[0]) + 1;
            std::array<std::vector<int>, 4> slab;
            std::array<std::vector<uint64_t>, 4> prefix;
            for (int channel = 0; channel < 4; ++channel) {
                slab[channel].resize(size_t(image[0]) * image[1], 0);
                prefix[channel].resize(pitch * (size_t(image[1]) + 1), 0);
            }
            const auto interface_face = [&](int axis, int x, int y, int z) {
                Shape next{{x, y, z}};
                if (++next[axis] >= image[axis]) return false;
                if (synthetic_layers)
                    return is_pore(x, y, z) != is_pore(next[0], next[1], next[2]);
                const double a = model.Distance(x + 1, y + 1, z + 1);
                const double b = model.Distance(next[0] + 1, next[1] + 1, next[2] + 1);
                // Match statistical analysis: zero distance is not an interface.
                return (a < 0.0 && b > 0.0) || (a > 0.0 && b < 0.0);
            };
            auto add_plane = [&](int z, int sign) {
                for (int y = 0; y < image[1]; ++y)
                    for (int x = 0; x < image[0]; ++x) {
                        const size_t index = size_t(y) * image[0] + x;
                        slab[0][index] += sign * is_pore(x, y, z);
                        slab[1][index] += sign * interface_face(0, x, y, z);
                        slab[2][index] += sign * interface_face(1, x, y, z);
                    }
            };
            auto add_z_faces = [&](int z, int sign) {
                for (int y = 0; y < image[1]; ++y)
                    for (int x = 0; x < image[0]; ++x)
                        slab[3][size_t(y) * image[0] + x] += sign * interface_face(2, x, y, z);
            };
            for (int z = 0; z < size[2]; ++z) add_plane(z, 1);
            for (int z = 0; z < size[2] - 1; ++z) add_z_faces(z, 1);
            const auto rectangle_sum = [&](int channel, int x, int y, int width, int height) {
                const auto &table = prefix[channel];
                return (table[size_t(y + height) * pitch + x + width] -
                        table[size_t(y) * pitch + x + width]) -
                       (table[size_t(y + height) * pitch + x] -
                        table[size_t(y) * pitch + x]);
            };
            for (int z = 0; z <= image[2] - size[2]; ++z) {
                if (z > 0) {
                    add_plane(z - 1, -1);
                    add_plane(z + size[2] - 1, 1);
                    if (size[2] > 1) {
                        add_z_faces(z - 1, -1);
                        add_z_faces(z + size[2] - 2, 1);
                    }
                }
                for (int channel = 0; channel < 4; ++channel) {
                    for (int y = 0; y < image[1]; ++y) {
                        uint64_t row = 0;
                        for (int x = 0; x < image[0]; ++x) {
                            row += slab[channel][size_t(y) * image[0] + x];
                            prefix[channel][size_t(y + 1) * pitch + x + 1] =
                                prefix[channel][size_t(y) * pitch + x + 1] + row;
                        }
                    }
                }
                for (int y = 0; y <= image[1] - size[1]; ++y) {
                    for (int x = 0; x <= image[0] - size[0]; ++x) {
                        const uint64_t count = rectangle_sum(0, x, y, size[0], size[1]);
                        const uint64_t faces = rectangle_sum(1, x, y, size[0] - 1, size[1]) +
                                               rectangle_sum(2, x, y, size[0], size[1] - 1) +
                                               rectangle_sum(3, x, y, size[0], size[1]);
                        best.considerGeometry({{x, y, z}}, count / volume, stat.mean,
                                              (faces / volume) / h, stat.surface_mean);
                    }
                }
            }
        } else {
            int bc = domain->getWithDefault<int>("BC", 0);
            if (selection_flow_params) {
                if (selection_flow_params->keyExists("BC"))
                    bc = selection_flow_params->getScalar<int>("BC");
                else if (selection_flow_params->keyExists("BoundaryCondition"))
                    bc = selection_flow_params->getScalar<int>("BoundaryCondition");
            }
            if (bc != 0)
                throw std::runtime_error("Permeability selection requires periodic body-force flow (BC = 0)");
            const int max_steps = selection_flow_params
                ? selection_flow_params->getWithDefault<int>("timestepMax", 100000) : 100000;
            const int interval = selection_flow_params
                ? selection_flow_params->getWithDefault<int>("analysis_interval", 1000) : 1000;
            const double tolerance = selection_flow_params
                ? selection_flow_params->getWithDefault<double>("tolerance", 1e-8) : 1e-8;
            if (!device_ready) {
                ScaLBL_SetDevice(rank);
                device_ready = true;
            }
            Array<signed char> window(size[0], size[1], size[2]);
            for (int z = 0; z <= image[2] - size[2]; ++z) {
                for (int y = 0; y <= image[1] - size[1]; ++y) {
                    for (int x = 0; x <= image[0] - size[0]; ++x) {
                        for (int k = 0; k < size[2]; ++k)
                            for (int j = 0; j < size[1]; ++j)
                                for (int i = 0; i < size[0]; ++i)
                                    window(i, j, k) = is_pore(x + i, y + j, z + k) ? 1 : 0;
                        rev_selection::HydraulicProperties properties;
                        try {
                            properties = rev_selection::simulateHydraulicProperties(
                                window, h, model.tau, {{model.Fx, model.Fy, model.Fz}},
                                max_steps, interval, tolerance, comm);
                        } catch (const std::exception &error) {
                            std::ostringstream message;
                            message << "Permeability candidate (" << x << ',' << y << ',' << z
                                    << "): " << error.what();
                            throw std::runtime_error(message.str());
                        }
                        best.considerHydraulics({{x, y, z}}, properties.permeability, stat.mean,
                                                properties.tortuosity, stat.tortuosity_mean);
                    }
                }
                std::cout << "Permeability selection: " << best.candidates << '/' << expected
                          << " windows evaluated" << std::endl;
            }
        }
        if (best.candidates != expected)
            throw std::runtime_error("Subdomain selection did not evaluate every valid window");
        if (!best.found) {
            write_status(name, "no_valid_tortuosity", &stat, best.candidates);
            std::cout << "Hydraulic selection found no candidate with defined tortuosity after "
                      << best.candidates << " windows." << std::endl;
            continue;
        }
        const std::string &raw_file = raw_files[property];
        if (!raw_file.empty())
            source->save(raw_file, best.start, size);
        const std::string domain_file = porosity ? "geometrical-rev-domain.db" : "hydraulic-rev-domain.db";
        source->saveDomain(domain_file, *domain, best.start, size, h);
        const double error = best.deviation;
        summary << name << ",selected," << stat.size << ',' << stat.size * h;
        for (int d : size) summary << ',' << d;
        for (int d : size) summary << ',' << d * h;
        for (int d : best.start) summary << ',' << d;
        for (int d = 0; d < 3; ++d) summary << ',' << best.start[d] + size[d];
        summary << ',' << stat.mean << ',' << best.value << ',' << error << ',' << best.candidates
                << ',' << csvString(raw_file) << ',' << csvString(source->filename)
                << ',' << csvString(source->read_type) << ',' << source->bytes;
        for (int d : source->shape) summary << ',' << d;
        for (int d = 0; d < 3; ++d) summary << ',' << best.start[d] + source->offset[d];
        for (int d = 0; d < 3; ++d) summary << ',' << best.start[d] + source->offset[d] + size[d];
        summary << ',' << csvString(domain_file);
        if (porosity)
            summary << ',' << stat.surface_mean << ',' << best.secondary_value << ',' << best.secondary_deviation
                    << ',' << best.score << ',' << porosity_weight << ',' << surface_weight << ",,,,,";
        else
            summary << ",,,," << best.score << ",,," << stat.tortuosity_mean << ',' << best.secondary_value
                    << ',' << best.secondary_deviation << ',' << permeability_weight << ',' << tortuosity_weight;
        summary << '\n';
        summary.flush();
        if (!summary) throw std::runtime_error("Failed to write subdomain_selection.csv");
        std::cout << "Selected " << name << " window at (" << best.start[0] << ','
                  << best.start[1] << ',' << best.start[2] << "), value=" << best.value
                  << ", relative deviation=" << error;
        if (porosity)
            std::cout << ", specific surface=" << best.secondary_value;
        else
            std::cout << ", tortuosity=" << best.secondary_value;
        std::cout << ", selection score=" << best.score;
        std::cout << std::endl;
    }
}

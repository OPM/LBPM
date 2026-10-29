#include "analysis/rev_subdomain_flow.h"

#include "common/Database.h"
#include "common/Domain.h"
#include "common/ScaLBL.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

namespace rev_selection {
namespace {

// Simulation storage must be released after every window, including failures.
template <class T> class DeviceBuffer {
public:
    explicit DeviceBuffer(size_t count) {
        void *storage = nullptr;
        ScaLBL_AllocateDeviceMemory(&storage, count * sizeof(T));
        data_ = static_cast<T *>(storage);
        if (!data_)
            throw std::runtime_error("Cannot allocate REV window flow storage");
    }
    ~DeviceBuffer() { ScaLBL_FreeDeviceMemory(data_); }
    DeviceBuffer(const DeviceBuffer &) = delete;
    DeviceBuffer &operator=(const DeviceBuffer &) = delete;
    T *data() const { return data_; }

private:
    T *data_ = nullptr;
};

} // namespace

HydraulicProperties simulateHydraulicProperties(
    const Array<signed char> &labels, double voxel_length, double tau,
    std::array<double, 3> force, int max_steps, int interval, double tolerance,
    const Utilities::MPI &comm) {
    if (comm.getSize() != 1)
        throw std::runtime_error("REV window flow requires one MPI process");
    if (!std::isfinite(voxel_length) || voxel_length <= 0.0 ||
        !std::isfinite(tau) || tau <= 0.5 ||
        !std::isfinite(tolerance) || tolerance <= 0.0 || tolerance >= 1.0 ||
        max_steps < 2 || interval < 1)
        throw std::runtime_error("Invalid REV window flow simulation settings");
    for (double component : force)
        if (!std::isfinite(component))
            throw std::runtime_error("REV window flow force must be finite");
    const double force_magnitude =
        std::hypot(std::hypot(force[0], force[1]), force[2]);
    if (!std::isfinite(force_magnitude) || force_magnitude <= 0.0)
        throw std::runtime_error("REV window flow requires a nonzero body force");

    const size_t limit = size_t(std::numeric_limits<int>::max());
    size_t padded_volume = 1;
    std::array<int, 3> size;
    for (int axis = 0; axis < 3; ++axis) {
        const size_t length = labels.size(axis);
        if (length == 0 || length > limit - 2 ||
            padded_volume > limit / (length + 2))
            throw std::runtime_error("REV window exceeds flow solver size limits");
        padded_volume *= length + 2;
        size[axis] = static_cast<int>(length);
    }
    if (labels.ndim() > 3 || labels.length() > limit / 19 - 32)
        throw std::runtime_error("REV window exceeds D3Q19 indexing limits");

    int pore_count = 0;
    for (size_t n = 0; n < labels.length(); ++n)
        pore_count += labels(n) > 0 ? 1 : 0;
    if (pore_count == 0)
        return {0.0, std::numeric_limits<double>::quiet_NaN()};

    auto domain_db = std::make_shared<Database>();
    domain_db->putVector<int>("n", {size[0], size[1], size[2]});
    domain_db->putVector<int>("nproc", {1, 1, 1});
    domain_db->putScalar<int>("BC", 0);
    domain_db->putScalar<double>("voxel_length", voxel_length);
    auto domain = std::make_shared<Domain>(domain_db, comm);

    // Match periodic MRT boundaries, including the geometry halo. No labels in
    // the source window are changed and no inlet/outlet reservoirs are added.
    const int nx = domain->Nx, ny = domain->Ny, nz = domain->Nz;
    for (int k = 0; k < nz; ++k) {
        const int z = (k + size[2] - 1) % size[2];
        for (int j = 0; j < ny; ++j) {
            const int y = (j + size[1] - 1) % size[1];
            for (int i = 0; i < nx; ++i) {
                const int x = (i + size[0] - 1) % size[0];
                domain->id[(k * ny + j) * nx + i] = labels(x, y, z);
            }
        }
    }
    domain->CommInit();
    auto communicator = std::make_shared<ScaLBL_Communicator>(domain);
    IntArray map(nx, ny, nz);
    const int allocation_count = (pore_count / 16 + 2) * 16;
    std::vector<int> neighbors(size_t(18) * allocation_count, 0);
    const int np = communicator->MemoryOptimizedLayoutAA(
        map, neighbors.data(), domain->id.data(), pore_count, 1);
    if (np > allocation_count)
        throw std::runtime_error("REV window flow layout exceeds allocation");

    DeviceBuffer<int> device_neighbors(size_t(18) * np);
    DeviceBuffer<double> distributions(size_t(19) * np);
    DeviceBuffer<double> velocity(size_t(3) * np);
    ScaLBL_CopyToDevice(device_neighbors.data(), neighbors.data(),
                        size_t(18) * np * sizeof(int));
    ScaLBL_D3Q19_Init(distributions.data(), np);
    const double relaxation_a = 1.0 / tau;
    const double relaxation_b =
        8.0 * (2.0 - relaxation_a) / (8.0 - relaxation_a);
    // The momentum estimator reads post-collision j/rho - F/2. Initialize
    // j=F/2 so the physical initial velocity is zero, including sealed pores
    // whose bounce-back otherwise preserves an artificial odd/even oscillation.
    ScaLBL_D3Q19_AAeven_MRT(
        distributions.data(), 0, np, np, relaxation_a, relaxation_b,
        0.5 * force[0], 0.5 * force[1], 0.5 * force[2]);
    ScaLBL_DeviceBarrier();

    std::vector<double> host_velocity(size_t(3) * np);
    const double viscosity = (tau - 0.5) / 3.0;
    const double porosity = double(pore_count) / double(labels.length());
    bool have_previous = false;
    double previous_flow = 0.0, previous_speed = 0.0, previous_z = 0.0;
    double previous_tortuosity = std::numeric_limits<double>::quiet_NaN();
    int stable_checks = 0, last_check = 0;

    // Use the same odd/even AA updates and second-order momentum estimator as
    // ScaLBL_MRTModel::Run. Check only after complete pairs of timesteps.
    const int final_step = max_steps - max_steps % 2;
    for (int step = 0; step < final_step;) {
        communicator->SendD3Q19AA(distributions.data());
        ScaLBL_D3Q19_AAodd_MRT(
            device_neighbors.data(), distributions.data(),
            communicator->FirstInterior(), communicator->LastInterior(), np,
            relaxation_a, relaxation_b, force[0], force[1], force[2]);
        communicator->RecvD3Q19AA(distributions.data());
        ScaLBL_D3Q19_AAodd_MRT(
            device_neighbors.data(), distributions.data(), 0,
            communicator->LastExterior(), np, relaxation_a, relaxation_b,
            force[0], force[1], force[2]);
        ScaLBL_DeviceBarrier();

        communicator->SendD3Q19AA(distributions.data());
        ScaLBL_D3Q19_AAeven_MRT(
            distributions.data(), communicator->FirstInterior(),
            communicator->LastInterior(), np, relaxation_a, relaxation_b,
            force[0], force[1], force[2]);
        communicator->RecvD3Q19AA(distributions.data());
        ScaLBL_D3Q19_AAeven_MRT(
            distributions.data(), 0, communicator->LastExterior(), np,
            relaxation_a, relaxation_b, force[0], force[1], force[2]);
        ScaLBL_DeviceBarrier();
        step += 2;
        if (step - last_check < interval && step != final_step)
            continue;
        last_check = step;

        ScaLBL_D3Q19_Momentum_2nd_order(
            distributions.data(), velocity.data(), np,
            force[0], force[1], force[2]);
        ScaLBL_DeviceBarrier();
        ScaLBL_CopyToHost(host_velocity.data(), velocity.data(),
                          size_t(3) * np * sizeof(double));
        double projected_sum = 0.0, speed_sum = 0.0, z_sum = 0.0;
        for (int k = 1; k < nz - 1; ++k)
            for (int j = 1; j < ny - 1; ++j)
                for (int i = 1; i < nx - 1; ++i) {
                    const int n = map(i, j, k);
                    if (n < 0)
                        continue;
                    std::array<double, 3> local_velocity;
                    for (int axis = 0; axis < 3; ++axis) {
                        const double component = host_velocity[size_t(axis) * np + n];
                        if (!std::isfinite(component))
                            throw std::runtime_error("Nonfinite REV window flow velocity");
                        projected_sum += component * force[axis] / force_magnitude;
                        local_velocity[axis] = component;
                    }
                    speed_sum += std::hypot(std::hypot(local_velocity[0], local_velocity[1]),
                                             local_velocity[2]);
                    z_sum += local_velocity[2];
                }
        const double flow = projected_sum / pore_count;
        const double mean_speed = speed_sum / pore_count;
        const double mean_z = z_sum / pore_count;
        if (!std::isfinite(flow) || !std::isfinite(mean_speed) || !std::isfinite(mean_z))
            throw std::runtime_error("Nonfinite REV window flow measurement");
        // Recover the same dimensionless tortuosity as statistical REV from
        // this velocity field. Exclude unresolved z-flow residuals rather than
        // interpreting roundoff in sealed or transverse-flow windows as flow.
        const double zero_flow = tolerance * std::max(force_magnitude, mean_speed);
        double tortuosity = std::numeric_limits<double>::quiet_NaN();
        if (std::abs(mean_z) > zero_flow && mean_speed > zero_flow) {
            tortuosity = mean_speed / mean_z;
            if (!std::isfinite(tortuosity))
                throw std::runtime_error("Nonfinite REV window tortuosity");
        }

        // A force-sized floor permits sealed, zero-flow pore space to converge
        // without division by zero. Require two consecutive stable intervals.
        const auto stable = [&](double value, double previous) {
            const double scale = std::max(force_magnitude,
                                         std::max(std::abs(value), std::abs(previous)));
            return std::abs(value - previous) <= tolerance * scale;
        };
        const bool stable_tortuosity =
            (std::isnan(tortuosity) && std::isnan(previous_tortuosity)) ||
            (std::isfinite(tortuosity) && std::isfinite(previous_tortuosity) &&
             std::abs(tortuosity - previous_tortuosity) <= tolerance *
                 std::max(1.0, std::max(std::abs(tortuosity), std::abs(previous_tortuosity))));
        if (have_previous && stable(flow, previous_flow) &&
            stable(mean_speed, previous_speed) && stable(mean_z, previous_z) && stable_tortuosity)
            ++stable_checks;
        else
            stable_checks = 0;
        if (stable_checks >= 2) {
            if (flow < -tolerance * force_magnitude)
                throw std::runtime_error("REV window flow opposes its driving force");
            const double permeability = 1013.0 * voxel_length * voxel_length *
                                        viscosity * porosity *
                                        std::max(0.0, flow) / force_magnitude;
            if (!std::isfinite(permeability))
                throw std::runtime_error("Nonfinite REV window permeability");
            return {permeability, tortuosity};
        }
        previous_flow = flow;
        previous_speed = mean_speed;
        previous_z = mean_z;
        previous_tortuosity = tortuosity;
        have_previous = true;
    }
    throw std::runtime_error("REV window permeability simulation did not converge within timestepMax");
}

double simulatePermeability(const Array<signed char> &labels,
                            double voxel_length, double tau,
                            std::array<double, 3> force, int max_steps,
                            int interval, double tolerance,
                            const Utilities::MPI &comm) {
    return simulateHydraulicProperties(labels, voxel_length, tau, force,
                                        max_steps, interval, tolerance, comm).permeability;
}

} // namespace rev_selection

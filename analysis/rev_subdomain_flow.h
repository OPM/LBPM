#ifndef LBPM_REV_SUBDOMAIN_FLOW_H
#define LBPM_REV_SUBDOMAIN_FLOW_H

#include "common/Array.h"
#include "common/MPI.h"

#include <array>

namespace rev_selection {

struct HydraulicProperties {
    double permeability;
    double tortuosity;
};

// Run a fresh, periodic, body-force-driven D3Q19 MRT solve on a window.
// Labels contain the interior only, in x-fast order; positive labels are pores.
// The returned permeability uses the same mDarcy conversion as ScaLBL_MRTModel.
// Tortuosity matches the displayed statistical REV quantity sum(|u|)/sum(u_z).
// It is NaN when the window has no resolved net z flow, including solid and
// sealed pore space. Its sign follows the existing z-velocity convention.
// A nonfinite or unconverged simulation raises std::runtime_error.
HydraulicProperties simulateHydraulicProperties(
    const Array<signed char> &labels, double voxel_length, double tau,
    std::array<double, 3> force, int max_steps, int interval, double tolerance,
    const Utilities::MPI &comm);

// Compatibility interface for callers that only need permeability.
double simulatePermeability(const Array<signed char> &labels,
                            double voxel_length, double tau,
                            std::array<double, 3> force, int max_steps,
                            int interval, double tolerance,
                            const Utilities::MPI &comm);

} // namespace rev_selection

#endif

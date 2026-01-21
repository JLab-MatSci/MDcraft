#pragma once

#include <mdcraft/configuration.h>

#include <vector>
#include <mdcraft/data/vector.h>
#include <mdcraft/data/matrix.h>

#ifdef DEBUG
#include <mdcraft/tools/helpers.h>
using namespace mdcraft::tools;
#endif

#ifdef mdcraft_ENABLE_BPNN
#include <unsupported/Eigen/CXX11/Tensor> 
#endif

#ifdef mdcraft_ENABLE_CUDA_THRUST
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#endif

#include <iostream>

namespace mdcraft {
namespace data {

/// @brief Simple structure for representing atom properties.
/// The structure is 8-byte aligned (take this into account when adding new fields).
struct Atom {
    vector           r;    // coordinate
    vector           v;    // velocity
    vector           f;    // force
    matrix           V;    // virial
    double           m;    // mass
    double           rcut; // radius of interaction
    double           rns;  // radius for neibs search
    double           Ep;   // potential energy
    double           Em;   // embedding energy (for EAM)
    double           nc;   // charge density (for EAM)
    unsigned short   bc;   // boundary condition type
    unsigned short kind;   // atom type
    unsigned int    uid;   // atom unique id

    // atom tag (zero value is reserved to early return from functions)
    unsigned short   tag = std::numeric_limits<unsigned short>::max();
};

using Atoms = std::vector<Atom>;

#ifdef mdcraft_ENABLE_CUDA_THRUST
using AtomsDev = thrust::device_vector<Atom>;
#endif

} // data
} // mdcraft

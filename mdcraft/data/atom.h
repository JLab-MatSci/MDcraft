#pragma once
#include <vector>
#include <functional>
#include <limits>
#include <mdcraft/data/vector.h>
#include <mdcraft/data/matrix.h>

#ifdef DEBUG
#include <mdcraft/tools/helpers.h>
using namespace mdcraft::tools;
#endif

#ifdef mdcraft_ENABLE_BPNN
#include <unsupported/Eigen/CXX11/Tensor> 
#endif

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
    double           nc_mod[2];   // charge density (for EAM)
    unsigned short   bc;   // boundary condition type
    unsigned short kind;   // atom type
    unsigned int    uid;   // atom unique id

    // atom tag (zero value is reserved to early return from functions)
    unsigned short   tag = std::numeric_limits<unsigned short>::max();
};

using Atoms = std::vector<Atom>;

// Constant-time ownership check, including when a pass contains no neighbours.
inline bool owns_atom(const Atoms& atoms, const Atom& atom) {
    const std::less<const Atom*> before;
    return !atoms.empty() && !before(&atom, atoms.data()) &&
           before(&atom, atoms.data() + atoms.size());
}


} // data
} // mdcraft

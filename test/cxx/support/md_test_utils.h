#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <mdcraft/data/atom.h>
#include <mdcraft/lattice/domain.h>
#include <mdcraft/neibs/verlet-list.h>
#include <mdcraft/solver/single.h>
#include <mdcraft/solver/stepper/verlet.h>
#include <mdcraft/tools/filesystem.h>

namespace mdcraft::test {

using data::Atom;
using data::Atoms;
using data::matrix;
using data::vector;
using lattice::Domain;
using neibs::VerletList;
using solver::Single;
using solver::stepper::Verlet;
using tools::Threads;

inline constexpr double kAngstromToNm = 0.1;
inline constexpr double kBoltzmann = 0.008314462175;

struct LammpsData {
  Domain domain;
  Atoms atoms;
};

struct ThermodynamicState {
  double potential_energy = 0.0;
  double embedding_energy = 0.0;
  double kinetic_energy = 0.0;
  double total_energy = 0.0;
  double virial = 0.0;
  vector com_velocity = vector::Zero();
};

inline tools::filesystem::path source_root() {
  return tools::filesystem::path(MDCRAFT_SOURCE_DIR);
}

inline tools::filesystem::path argon_ball_path() {
  return source_root() / "test" / "cxx" / "LMP_BALL_130k_Ar";
}

inline tools::filesystem::path aluminum_eam_path() {
  return source_root() / "test" / "python" / "solver" / "Al2009.eam.alloy";
}

inline vector make_vector(double x, double y, double z) {
  vector v;
  v << x, y, z;
  return v;
}

inline LammpsData load_lammps_atomic_data(const tools::filesystem::path& filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open LAMMPS data file: " + filename.string());
  }

  std::string line;
  std::size_t total_atoms = 0;
  double xmin = 0.0;
  double xmax = 0.0;
  double ymin = 0.0;
  double ymax = 0.0;
  double zmin = 0.0;
  double zmax = 0.0;
  std::unordered_map<unsigned short, double> mass_by_type;

  while (std::getline(file, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }

    if (line.find("atoms") != std::string::npos && line.find("atom") != std::string::npos) {
      std::istringstream iss(line);
      iss >> total_atoms;
      continue;
    }

    if (line.find("xlo xhi") != std::string::npos) {
      std::istringstream iss(line);
      iss >> xmin >> xmax;
      continue;
    }
    if (line.find("ylo yhi") != std::string::npos) {
      std::istringstream iss(line);
      iss >> ymin >> ymax;
      continue;
    }
    if (line.find("zlo zhi") != std::string::npos) {
      std::istringstream iss(line);
      iss >> zmin >> zmax;
      continue;
    }

    if (line.find("Masses") != std::string::npos) {
      while (std::getline(file, line) && (line.empty() || line[0] == '#')) {
      }

      while (!line.empty() && line[0] != '#') {
        std::istringstream iss(line);
        unsigned short type = 0;
        double mass = 0.0;
        if (iss >> type >> mass) {
          mass_by_type[type] = mass;
        }
        if (!std::getline(file, line)) {
          break;
        }
      }
      continue;
    }

    if (line.find("Atoms") != std::string::npos) {
      while (std::getline(file, line) && (line.empty() || line[0] == '#')) {
      }
      break;
    }
  }

  Atoms atoms;
  atoms.reserve(total_atoms);

  do {
    if (line.empty() || line[0] == '#') {
      continue;
    }

    std::istringstream iss(line);
    Atom atom{};
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    if (!(iss >> atom.uid >> atom.kind >> x >> y >> z)) {
      continue;
    }

    atom.r << x * kAngstromToNm, y * kAngstromToNm, z * kAngstromToNm;
    const auto it = mass_by_type.find(atom.kind);
    atom.m = it == mass_by_type.end() ? 0.0 : it->second;
    atoms.push_back(atom);
  } while (std::getline(file, line) && atoms.size() < total_atoms);

  Domain domain(
    xmin * kAngstromToNm, xmax * kAngstromToNm,
    ymin * kAngstromToNm, ymax * kAngstromToNm,
    zmin * kAngstromToNm, zmax * kAngstromToNm
  );

  return {domain, std::move(atoms)};
}

inline void set_common_atom_properties(
  Atoms& atoms,
  double mass,
  double rcut,
  double neighbor_skin_factor
) {
  for (std::size_t i = 0; i < atoms.size(); ++i) {
    atoms[i].uid = atoms[i].uid == 0 ? static_cast<unsigned int>(i + 1) : atoms[i].uid;
    atoms[i].kind = atoms[i].kind == 0 ? static_cast<unsigned short>(1) : atoms[i].kind;
    atoms[i].m = mass;
    atoms[i].rcut = rcut;
    atoms[i].rns = neighbor_skin_factor * rcut;
  }
}

inline void assign_deterministic_velocities(Atoms& atoms, double scale) {
  vector momentum = vector::Zero();
  double total_mass = 0.0;

  for (std::size_t i = 0; i < atoms.size(); ++i) {
    const double sx = static_cast<double>(static_cast<int>(i % 7) - 3);
    const double sy = static_cast<double>(static_cast<int>((i / 7) % 5) - 2);
    const double sz = static_cast<double>(static_cast<int>((i / 35) % 3) - 1);
    atoms[i].v << scale * sx, scale * sy, scale * sz;
    momentum += atoms[i].m * atoms[i].v;
    total_mass += atoms[i].m;
  }

  vector com_velocity = vector::Zero();
  if (total_mass > 0.0) {
    com_velocity = momentum / total_mass;
  }
  for (auto& atom : atoms) {
    atom.v -= com_velocity;
  }
}

inline void zero_velocities(Atoms& atoms) {
  for (auto& atom : atoms) {
    atom.v = vector::Zero();
  }
}

inline ThermodynamicState measure_state(const Atoms& atoms) {
  ThermodynamicState state;
  double total_mass = 0.0;

  for (const auto& atom : atoms) {
    state.potential_energy += atom.Ep;
    state.embedding_energy += atom.Em;
    state.virial += atom.V.trace();
    state.com_velocity += atom.m * atom.v;
    total_mass += atom.m;
  }

  if (total_mass > 0.0) {
    state.com_velocity /= total_mass;
  } else {
    state.com_velocity = vector::Zero();
  }

  for (const auto& atom : atoms) {
    const vector relative_velocity = atom.v - state.com_velocity;
    state.kinetic_energy += 0.5 * atom.m * relative_velocity.squaredNorm();
  }

  state.total_energy = state.potential_energy + state.embedding_energy + state.kinetic_energy;
  return state;
}

inline ThermodynamicState recompute_state(
  Single& solver,
  Atoms& atoms,
  VerletList& nlist
) {
  auto& raw_list = nlist.get();
  solver.prepare(atoms, atoms, raw_list);
  for (auto& stage : solver.stages) {
    stage(atoms, atoms, raw_list);
  }
  solver.virials(atoms, atoms, raw_list);
  return measure_state(atoms);
}

inline std::vector<ThermodynamicState> run_md_steps(
  Single& solver,
  Verlet& stepper,
  Atoms& atoms,
  VerletList& nlist,
  double dt,
  int steps,
  int neighbor_update_interval
) {
  std::vector<ThermodynamicState> history;
  history.reserve(static_cast<std::size_t>(steps + 1));

  nlist.update(-1.0);
  history.push_back(recompute_state(solver, atoms, nlist));

  for (int step = 1; step <= steps; ++step) {
    stepper.make_step(atoms, nlist.get(), dt);
    if (neighbor_update_interval > 0 && (step % neighbor_update_interval) == 0) {
      nlist.update(-1.0);
    }
    history.push_back(recompute_state(solver, atoms, nlist));
  }

  return history;
}

inline std::size_t find_atom_index_by_uid(const Atoms& atoms, unsigned int uid) {
  const auto it = std::find_if(atoms.begin(), atoms.end(), [uid](const Atom& atom) {
    return atom.uid == uid;
  });
  if (it == atoms.end()) {
    throw std::runtime_error("Atom uid not found: " + std::to_string(uid));
  }
  return static_cast<std::size_t>(std::distance(atoms.begin(), it));
}

inline std::pair<Atoms, Domain> make_fcc_bulk(
  int nx,
  int ny,
  int nz,
  double lattice_constant,
  double mass
) {
  Atoms atoms;
  atoms.reserve(static_cast<std::size_t>(4 * nx * ny * nz));

  const vector basis[4] = {
    make_vector(0.0, 0.0, 0.0),
    make_vector(0.5, 0.5, 0.0),
    make_vector(0.5, 0.0, 0.5),
    make_vector(0.0, 0.5, 0.5),
  };

  const double lx = nx * lattice_constant;
  const double ly = ny * lattice_constant;
  const double lz = nz * lattice_constant;
  const vector origin = make_vector(-0.5 * lx, -0.5 * ly, -0.5 * lz);

  unsigned int uid = 1;
  for (int iz = 0; iz < nz; ++iz) {
    for (int iy = 0; iy < ny; ++iy) {
      for (int ix = 0; ix < nx; ++ix) {
        const vector cell_origin = make_vector(
          origin.x() + ix * lattice_constant,
          origin.y() + iy * lattice_constant,
          origin.z() + iz * lattice_constant
        );

        for (const auto& shift : basis) {
          Atom atom{};
          atom.uid = uid++;
          atom.kind = 1;
          atom.m = mass;
          atom.r = cell_origin + lattice_constant * shift;
          atoms.push_back(atom);
        }
      }
    }
  }

  Domain domain(
    -0.5 * lx, 0.5 * lx,
    -0.5 * ly, 0.5 * ly,
    -0.5 * lz, 0.5 * lz
  );
  domain.set_periodic(0, true);
  domain.set_periodic(1, true);
  domain.set_periodic(2, true);
  return {std::move(atoms), domain};
}

}  // namespace mdcraft::test

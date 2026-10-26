#include <cmath>
#include <mdcraft/solver/single.h>

namespace mdcraft {
namespace solver {

using ::mdcraft::tools::Threads;
using ::mdcraft::data::vector;
using ::mdcraft::data::matrix;

Single::Single(
	Potential&  potential,
	Boundary&   boundary,
	Thermostat& thermostat,
	Threads&    threads
) : Pair(
		boundary,
	    thermostat,
	    threads 
    ), 
    potential(potential)
{
	if (potential.type == potential::type::pair) {
		stages.emplace_back([this] (Atoms& atoms, Atoms& neibs, NeibsList& nlist) {
			forces(atoms, neibs, nlist);
		});
	}
	if (potential.type == potential::type::eam) {
		stages.emplace_back([this] (Atoms& atoms, Atoms& neibs, NeibsList& nlist) {
			charge(atoms, neibs, nlist);
		});
		stages.emplace_back([this] (Atoms& atoms, Atoms& neibs, NeibsList& nlist) {
			forces(atoms, neibs, nlist);
		});
	}
}

void Single::prepare_one(
	Atoms::iterator atomit, 
	Atoms&          neibs, 
	NeibsListOne&   nlist
) {
	auto& atom = *atomit;
	atom.f  = vector::Zero();
	atom.V  = matrix::Zero();
	atom.Ep = 0.0;
	if (potential.type == potential::type::eam){
		atom.Em = 0.0;
		atom.nc = 0.0;
	}
	for (auto& rho : atom.nc_mod) rho = 0.0;
	atom.rcut = potential.rcut();
}

void Single::charge_one(
	Atoms::iterator atomit,
	Atoms&          neibs, 
	NeibsListOne&   nlist
){
	auto& atom = *atomit;
	if (atom.tag == 0) return;
	potential.density(atomit, neibs, nlist);
}

void Single::force_one(
	Atoms::iterator atomit,
	Atoms&          neibs, 
	NeibsListOne&   nlist
) {
	auto& atom = *atomit;
	if (atom.tag == 0) return;
	potential.force(atomit, neibs, nlist);
}

void Single::virial_one(
	Atoms::iterator atomit,
	Atoms&          neibs, 
	NeibsListOne&   nlist
) {
	auto& atom = *atomit;
	if (atom.tag == 0) return;
	potential.virial(atomit, neibs, nlist);
}

} // solver
} // mdcraft

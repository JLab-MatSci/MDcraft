#pragma once

#include "pair.h"
#include <mdcraft/solver/potential/cuda/base_potential.cuh>

namespace mdcraft::solver::cuda {

class Single : public TPairSolver<Single>
{
public:
	Single() {}
	// Single(BasePotential& pot) {}

	void forces_(Atoms& atoms, Atoms& neibs, NeibsList& nlist)
	{
		static int i = 0;
		std::cout << ++i << std::endl;
		// forces(atoms, neibs, nlist);
		// using TPairSolver<Single>::forces;
		/*TPairSolver<Single>::*/
		// static_cast<> forces(atoms, neibs, nlist);
		//forces(atoms, neibs, nlist);
	}
};

template <typename P>
class TSingle : public Single
{
public:
	// Single(
	// 	Boundary&   boundary   = boundary::dummy_boundary,
	// 	Thermostat& thermostat = thermostat::dummy_thermostat
	// );
	TSingle() {}

	void forces(Atoms& atoms, Atoms& neibs, NeibsList& nlist) /*const*/
	{
		mdcraft::solver::potential::cuda::k_forces<P><<<1,1>>>(potential_.get());
	}

private:
	std::shared_ptr<Potential<P>> potential_;
};

} // namespace mdcraft::solver::cuda

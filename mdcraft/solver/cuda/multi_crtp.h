#pragma once

#include "pair.h"
#include <mdcraft/solver/potential/cuda/base_potential.cuh>

#include <mdcraft/solver/thermostat/base.h>
#include <mdcraft/solver/boundary/base.h>

namespace mdcraft::solver::cuda {

class Multi : public TPairSolver<Multi>
{
public:
	Multi() {}
	Multi(BasePotential& pot) {}
};

template <typename P>
class TMulti : public Multi
{
public:
	TMulti() {}

	void forces(Atoms& atoms, Atoms& neibs, NeibsList& nlist) /*const*/
	{
		potential_.force();
	}

private:
	std::shared_ptr<Potential<P>> potential_;
};

} // namespace mdcraft::solver::cuda

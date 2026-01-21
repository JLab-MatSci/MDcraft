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
	}
};

template <typename P>
class TSingle : public Single
{
public:
	TSingle() {}

	void forces(Atoms& atoms, Atoms& neibs, NeibsList& nlist) /*const*/
	{
		mdcraft::solver::potential::cuda::k_forces<P><<<1,1>>>(potential_.get());
	}

private:
	std::shared_ptr<Potential<P>> potential_;
};

} // namespace mdcraft::solver::cuda

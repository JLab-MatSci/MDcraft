#pragma once

#include <mdcraft/configuration.h>

#include <mdcraft/data/atom.h>
#include <mdcraft/configuration.h>
#include <mdcraft/neibs/verlet-list.h>

#include <mdcraft/solver/potential/cuda/base_potential.cuh>

namespace mdcraft::solver::cuda {

using Atoms     = mdcraft::data::Atoms;
using AtomsDev  = mdcraft::data::AtomsDev;
using NeibsList = neibs::List;
using NeibsListOne = ::mdcraft::neibs::ListOne;

using BasePotential = mdcraft::solver::potential::cuda::BasePotential;

template <typename P>
void inline CUDA_ASSERT_POINTER_ON_DEVICE(P* p)
{
	cudaPointerAttributes attr;
	cudaError_t err = cudaPointerGetAttributes(&attr, p);
	if (attr.type != cudaMemoryTypeDevice) assert(false);
}

class BaseSolver
{
public:
	BaseSolver() = default;
	BaseSolver(BasePotential& pot) {}
};

template <typename S>
class TBaseSolver : public BaseSolver
{
public:
	TBaseSolver() = default;

	void prepare(Atoms& atoms, Atoms& neibs, NeibsList& nlist) /*const*/
	{
		static_cast<S*>(this)->prepare(atoms, neibs, nlist);
	}
};

} // namespace mdcraft::solver::cuda
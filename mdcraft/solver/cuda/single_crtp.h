#pragma once

#include "pair_crtp.h"
#include <mdcraft/solver/potential/cuda/base_potential.cuh>

namespace mdcraft::solver::cuda {

template <typename P>
class TSingle : public TPairSolver<TSingle<P>>
{
public:
	TSingle() = default;

	// for now ctor copies potential from host to device
	// it is strange, as only potential itself must be
	// responsible for its resource management
	TSingle(Potential<P>* p)
	{
		static_assert(std::is_trivially_copyable_v<std::remove_reference_t<P>>);
		cudaMalloc(&potential_dev_, sizeof(P));
		cudaMemcpy(potential_dev_, p, sizeof(P), cudaMemcpyHostToDevice);
	}

	// NO DTOR TO FREE MEMORY. ELSE WE LOOSE TRIVIALLY COPYABLE CLASS
	/*~TSingle()
	{
		cudaFree(potential_dev_);
	}*/

#ifdef mdcraft_ENABLE_CUDA_THRUST
	__device__ void force_one(
		data::Atom&             atom,
		data::Atom*             neibs,
		neibs::NeibsListOneDev& nlist
	)
	{
		potential_dev_->force(atom, neibs, nlist);
	}
#endif

private:
	Potential<P>* potential_dev_;

};

} // namespace mdcraft::solver::cuda

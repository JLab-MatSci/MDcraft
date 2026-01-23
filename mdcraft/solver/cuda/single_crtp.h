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
	TSingle(P* p)
	{
		P* p_dev = p->ptr_dev();
		CUDA_ASSERT_POINTER_ON_DEVICE(p_dev);
		potential_dev_ = p_dev;
	}

#ifdef mdcraft_ENABLE_CUDA_THRUST
	__device__ void force_one(
		data::Atom&             atom,
		data::Atom*             neibs,
		neibs::NeibsListOneDev& nlist
	)
	{
		potential_dev_->force_impl(atom, neibs, nlist);
	}
#endif

private:
	P* potential_dev_ = nullptr;

};

} // namespace mdcraft::solver::cuda

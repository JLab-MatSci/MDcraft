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

#ifdef mdcraft_ENABLE_CUDA_THRUST
#define DISPATCH_SOLVER_FOR_CUDA(T, device_ptr)                  \
	T* obj = static_cast<T*>(this);                              \
                                                                 \
	/* convert into pointer accessible on device */              \
	thrust::device_vector<T> v_obj_dev(1, *obj);                 \
	T* device_ptr = thrust::raw_pointer_cast(v_obj_dev.data());  \
	CUDA_ASSERT_POINTER_ON_DEVICE(device_ptr);
#endif

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
};

} // namespace mdcraft::solver::cuda
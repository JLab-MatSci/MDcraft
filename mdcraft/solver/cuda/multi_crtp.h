#pragma once

#include "pair_crtp.h"
#include <mdcraft/solver/potential/cuda/base_potential.cuh>

#include <mdcraft/solver/thermostat/base.h>
#include <mdcraft/solver/boundary/base.h>

namespace mdcraft::solver::cuda {

template <typename P>
class TMulti : public TPairSolver<TMulti<P>, P>
{
public:
	using potential_t = P;

	TMulti() = default;

	__device__ void force_one_impl(
		data::Atom&             atom,
		data::Atom*             neibs,
		neibs::NeibsListOneDev& nlist
	)
	{
	}

	// TEMPORALLY !!! Will be changed to some sort of array of ptrs.
	// IT COMPILES! EVEN THOUNGH THIS ONE RETURNS int, AND 
	// THE ONE IN TSingle RETURNS potential_t*
	// __host__ __device__ int potential_impl() /*const*/
	// {
	// 	return 0;
	// }
};

} // namespace mdcraft::solver::cuda

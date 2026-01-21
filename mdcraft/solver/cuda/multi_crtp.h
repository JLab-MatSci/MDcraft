#pragma once

#include "pair_crtp.h"
#include <mdcraft/solver/potential/cuda/base_potential.cuh>

#include <mdcraft/solver/thermostat/base.h>
#include <mdcraft/solver/boundary/base.h>

namespace mdcraft::solver::cuda {

class Multi : public TPairSolver<Multi>
{
public:
	Multi() = default;
};

template <typename P>
class TMulti : public Multi
{
public:
	TMulti() = default;

	__device__ void force_one_impl(
		data::Atom&             atom,
		data::Atom*             neibs,
		neibs::NeibsListOneDev& nlist
	)
	{
	}
};

} // namespace mdcraft::solver::cuda

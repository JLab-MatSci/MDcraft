#pragma once

#include "isolver.h"

namespace mdcraft::solver::cuda {

template <typename P>
using Potential = mdcraft::solver::potential::cuda::TBasePotential<P>;

class PairSolver : public TBaseSolver<PairSolver>
{
public:
	PairSolver() {}

	void prepare(Atoms& atoms, Atoms& neibs, NeibsList& nlist) /*const*/
	{

	}
};

template <typename S>
class TPairSolver : public PairSolver /*public TBaseSolver<Pair<S>>*/
{
public:
	TPairSolver() {}

	void forces(Atoms& atoms, Atoms& neibs, NeibsList& nlist) /*const*/
	{
		std::cout << "TPairSolver::forces" << std::endl;
		static_cast<S*>(this)->forces(atoms, neibs, nlist);
	}
};

} // namespace mdcraft::solver::cuda
#include <mdcraft/data/atom.h>
#include <mdcraft/configuration.h>
#include <mdcraft/neibs/verlet-list.h>

#include <mdcraft/solver/potential/cuda/base_potential.cuh>

namespace mdcraft::solver::cuda {

using Atoms = mdcraft::data::Atoms;
using NeibsList = neibs::List;

using BasePotential = mdcraft::solver::potential::cuda::BasePotential;

class BaseSolver
{
public:
	BaseSolver() {}
	BaseSolver(BasePotential& pot) {}
};

template <typename S>
class TBaseSolver : public BaseSolver
{
public:
	TBaseSolver() {}

	void prepare(Atoms& atoms, Atoms& neibs, NeibsList& nlist) /*const*/
	{
		static_cast<S*>(this)->prepare(atoms, neibs, nlist);
	}
};

} // namespace mdcraft::solver::cuda
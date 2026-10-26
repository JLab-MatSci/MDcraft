#pragma once

#include <mdcraft/configuration.h>
#include <mdcraft/data/atom.h>
#include <mdcraft/neibs/verlet-list.h>

#ifdef mdcraft_ENABLE_MPI
#include <mdcraft/decomp/decomp.h>
#endif

namespace mdcraft::solver {

using data::Atoms;
using NeibsList = neibs::List;
#ifdef mdcraft_ENABLE_MPI
using decomp::Decomp;
#endif

using StageFunc = std::function<void(Atoms&, Atoms&, NeibsList&)>;

class ISolver {
public:
	virtual ~ISolver();
	virtual void evaluate(Atoms& atoms, NeibsList& nlist);
#ifdef mdcraft_ENABLE_MPI
	virtual void evaluate(Decomp& decomp, NeibsList& local, NeibsList& alien);
#endif

	virtual void prepare(Atoms& atoms, Atoms& neibs, NeibsList& nlist);

#ifdef mdcraft_ENABLE_MPI
	virtual void prepare(Decomp& decomp, NeibsList& nlist1, NeibsList& nlist2);
#endif

	std::vector<StageFunc> stages;
};

} // namespace mdcraft::solver
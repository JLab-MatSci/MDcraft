#include <mdcraft/solver/boundary/periodic.h>

namespace mdcraft {
namespace solver {
namespace boundary {

Periodic::Periodic(Domain& domain): Base(domain) {}

unsigned short Periodic::need_apply(
	Atoms::iterator atomit
) {
	auto& atom = *atomit;
	if (atom.tag == 0) return 0x0;

	bool need = false;
	for (auto d = 0; d < 3; ++d) {
		need = need || domain.nearest_distance(atom.r, d) < atom.rns;
	}

	return need ? type : 0x0;
}

void Periodic::prepare_one(
	Atoms::iterator atomit, 
	Atoms&          neibs,
	NeibsListOne&   nlist
) {
	auto& atom = *atomit;
	atom.bc |= need_apply(atomit);
	if (!(atom.bc & type)) return;
	// auto ghosts = prepare_ghosts(atomit, neibs, nlist);
	// if (ghosts.size() == 0) return;
	// solver->prepare_one(particle, ghosts);
}

Atoms Periodic::prepare_ghosts(
	Atoms::iterator atomit, 
	Atoms&          neibs,
	NeibsListOne&   nlist
) {
	return make_ghosts(atomit, neibs, nlist, true);
}

Atoms Periodic::make_ghosts(
	Atoms::iterator atomit, Atoms& neibs, NeibsListOne& nlist,
	bool with_density
) {
	auto& atom = *atomit;
	Atoms ghosts; ghosts.reserve(2*nlist.size()/3);
	for (auto& j : nlist) {
		auto& neib = neibs[j];

		vector dist = neib.r - atom.r;

		auto r2rns = std::max(neib.rns, atom.rns); r2rns *= r2rns;

		if (dist.dot(dist) < r2rns) continue;

		dist = domain.shortest(dist);

		if (dist.dot(dist) > r2rns) continue;

		// Other workers may be updating neib.f, V, Ep/Em, or (during
		// charge) its densities. A whole-Atom copy reads those fields and
		// races with their writers. Copy only immutable inputs of this stage.
		ghosts.emplace_back();
		auto& ghost = ghosts.back();
		ghost.r = atom.r + dist;
		ghost.v = neib.v;
		ghost.m = neib.m;
		ghost.rcut = neib.rcut;
		ghost.rns = neib.rns;
		ghost.bc = neib.bc;
		ghost.kind = neib.kind;
		ghost.uid = neib.uid;
		ghost.tag = neib.tag;
		ghost.f.setZero();
		ghost.V.setZero();
		ghost.Ep = ghost.Em = 0.0;
		ghost.nc = with_density ? neib.nc : 0.0;
		for (int k = 0; k < 2; ++k)
			ghost.nc_mod[k] = with_density ? neib.nc_mod[k] : 0.0;
	}
	return ghosts;
}

void Periodic::charge_one(
	Atoms::iterator atomit,
	Atoms&          neibs,
	NeibsListOne&   nlist
) {
	auto& atom = *atomit;
	if (!(atom.bc & type)) return;
	// Density kernels use geometry, not the unfinished neighbour densities.
	auto ghosts = make_ghosts(atomit, neibs, nlist, false);
	if (ghosts.size() == 0) return;
	NeibsListOne nlist_ghosts(ghosts.size()); 
	for (std::size_t i = 0; i < ghosts.size(); ++i) nlist_ghosts[i] = i;
	solver->charge_one(atomit, ghosts, nlist_ghosts);
}

void Periodic::force_one(
	Atoms::iterator atomit,
	Atoms&          neibs,
	NeibsListOne&   nlist
) {
	auto& atom = *atomit;
	if (!(atom.bc & type)) return;
	auto ghosts = prepare_ghosts(atomit, neibs, nlist);
	if (ghosts.size() == 0) return;
	NeibsListOne nlist_ghosts(ghosts.size());
	for (std::size_t i = 0; i < ghosts.size(); ++i) nlist_ghosts[i] = i;
	solver->force_one(atomit, ghosts, nlist_ghosts);
}

void Periodic::virial_one(
	Atoms::iterator atomit,
	Atoms&          neibs,
	NeibsListOne&     nlist
) {
	auto& atom = *atomit;
	if (!(atom.bc & type)) return;
	auto ghosts = prepare_ghosts(atomit, neibs, nlist);
	if (ghosts.size() == 0) return;
	NeibsListOne nlist_ghosts(ghosts.size()); 
	for (std::size_t i = 0; i < ghosts.size(); ++i) nlist_ghosts[i] = i;
	solver->virial_one(atomit, ghosts, nlist_ghosts);
}

} // boundary
} // solver
} // mdcraft

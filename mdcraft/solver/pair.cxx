#include <mdcraft/solver/pair.h>

namespace mdcraft::solver {
#ifdef mdcraft_ENABLE_MPI
void Pair::evaluate(Decomp& decomp, NeibsList& local, NeibsList& alien) {
    prepare(decomp, local, alien);
    if (stages.size() > 1) charge(decomp, local, alien);
    forces(decomp, local, alien);
}
#endif


Pair::Pair(
	Boundary&   boundary,
	Thermostat& thermostat,
	Threads&    threads
) : thermostat (thermostat),
    boundary   (boundary),
    pool       (threads)
{
	boundary.bind_solver(*this);
}

void Pair::prepare(
	Atoms&     atoms,
	Atoms&     neibs, 
	NeibsList& nlist
) {
	auto apply = [&](Atoms::iterator atomit) {
		auto& atom = *atomit;
		auto i = std::distance(atoms.begin(), atomit);
		prepare_one(atomit, neibs, nlist[i]);
		atom.bc = 0x0; boundary.prepare_one(atomit, neibs, nlist[i]);
 	};

	pool.for_each(atoms.begin(), atoms.end(), apply);
}

#ifdef mdcraft_ENABLE_MPI
void Pair::prepare(
	Decomp&    decomp,
	NeibsList& nlist1,
	NeibsList& nlist2
) {
    decomp.exchange_start();
    // Reset all per-step accumulators and bc only ONCE per local atom.
    prepare(decomp.locals(), decomp.locals(), nlist1);
    decomp.exchange_end();
    // Boundary preparation must also see alien neighbours, but calling
    // prepare() here would erase freshly initialized boundary flags again.
    auto& atoms = decomp.locals();
    auto& neibs = decomp.aliens();
    auto apply = [&](Atoms::iterator atomit) {
        auto i = std::distance(atoms.begin(), atomit);
        boundary.prepare_one(atomit, neibs, nlist2[i]);
    };
    pool.for_each(atoms.begin(), atoms.end(), apply);
}
#endif

void Pair::charge(
	Atoms&     atoms, 
	Atoms&     neibs, 
	NeibsList& nlist
) {
	auto apply = [&](Atoms::iterator atomit) {
		auto& atom = *atomit;
		auto i = std::distance(atoms.begin(), atomit);
		charge_one(atomit, neibs, nlist[i]);
		if (atom.bc) boundary.charge_one(atomit, neibs, nlist[i]);
 	};

	pool.for_each(atoms.begin(), atoms.end(), apply);
}

#ifdef mdcraft_ENABLE_MPI
void Pair::charge(
	Decomp&    decomp,
	NeibsList& nlist1,
	NeibsList& nlist2
	) {
	decomp.exchange_start();
	charge(decomp.locals(), decomp.locals(), nlist1);
	decomp.exchange_end();
	charge(decomp.locals(), decomp.aliens(), nlist2);
}
#endif

void Pair::accumulate_forces(
    Atoms& atoms,
    Atoms& neibs,
    NeibsList& nlist,
    bool apply_thermostat
) {
    auto apply = [&](Atoms::iterator atomit) {
        auto& atom = *atomit;
        auto i = std::distance(atoms.begin(), atomit);
        force_one(atomit, neibs, nlist[i]);
        if (atom.bc) boundary.force_one(atomit, neibs, nlist[i]);
        if (apply_thermostat && atom.tag) thermostat.apply_one(atomit);
    };
    pool.for_each(atoms.begin(), atoms.end(), apply);
}

void Pair::forces(Atoms& atoms, Atoms& neibs, NeibsList& nlist) {
    accumulate_forces(atoms, neibs, nlist, true);
}

#ifdef mdcraft_ENABLE_MPI
void Pair::forces(
	Decomp&    decomp,
	NeibsList& nlist1,
	NeibsList& nlist2
	) {
    decomp.exchange_start();
    accumulate_forces(decomp.locals(), decomp.locals(), nlist1, false);
    decomp.exchange_end();
    accumulate_forces(decomp.locals(), decomp.aliens(), nlist2, false);
    // The thermostat acts on the FULL force, once per integration step.
    auto& atoms = decomp.locals();
    pool.for_each(atoms.begin(), atoms.end(), [&](Atoms::iterator atomit) {
        if (atomit->tag) thermostat.apply_one(atomit);
    });
}
#endif

void Pair::virials(Atoms& atoms, Atoms& neibs, NeibsList& nlist) {
    accumulate_virials(atoms, neibs, nlist, true);
}

void Pair::accumulate_virials(
	Atoms&     atoms, 
	Atoms&     neibs, 
	NeibsList& nlist, bool reset
) {
	auto apply = [&](Atoms::iterator atomit) {
		auto& atom = *atomit;
        if (reset) { atom.V.setZero(); atom.Ep = 0.0; atom.Em = 0.0; }
		auto i = std::distance(atoms.begin(), atomit);
		virial_one(atomit, neibs, nlist[i]);
		if (atom.bc) boundary.virial_one(atomit, neibs, nlist[i]);
	};

	pool.for_each(atoms.begin(), atoms.end(), apply);
}

#ifdef mdcraft_ENABLE_MPI
void Pair::virials(
	Decomp&    decomp,
	NeibsList& nlist1,
	NeibsList& nlist2
	) {
	decomp.exchange_start();
	accumulate_virials(decomp.locals(), decomp.locals(), nlist1, true);
	decomp.exchange_end();
	accumulate_virials(decomp.locals(), decomp.aliens(), nlist2, false);
}
#endif

} // namespace mdcraft::solver

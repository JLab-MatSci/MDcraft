#include <cmath>
#include <mdcraft/configuration.h>
#include <mdcraft/tools/threads.h>
#include <mdcraft/solver/stepper/verlet.h>

namespace mdcraft::solver::stepper {

using data::vector;

Verlet::Verlet(
	ISolver& solver,
	Threads& threads
) : Base(
        solver, 
        threads
	) 
{}

void Verlet::step_coords(Atoms& atoms, double dt) {
	auto apply = [&](Atoms::iterator atomit) {
		step_one_coords(atomit, dt);
	};

	pool.for_each(atoms.begin(),  atoms.end(), apply);

}

void Verlet::step_velocities(Atoms& atoms, double dt) {
    auto apply = [&](Atoms::iterator atomit) {
        step_one_velocities(atomit, dt);
    };

    pool.for_each(atoms.begin(), atoms.end(), apply);
}

void Verlet::make_step(
	Atoms&     atoms, 
	NeibsList& nlist, 
	double     dt
) {
	step_coords(atoms, dt);
	solver.evaluate(atoms, nlist);
	step_velocities(atoms, dt);
}

#ifdef mdcraft_ENABLE_MPI
void Verlet::make_step(
	Decomp&     decomp,
	NeibsList&  nlist1,
	NeibsList&  nlist2,
	double      dt
) {
	step_coords(decomp.locals(), dt);
	solver.evaluate(decomp, nlist1, nlist2);
	step_velocities(decomp.locals(), dt);
}
#endif

void Verlet::step_one_coords(Atoms::iterator atomit, double dt) {
	auto& atom = *atomit;
	if (atom.tag == 0) return;
	vector a = atom.f; a /= atom.m;
	vector v = atom.v; v += a*dt*0.5;
	atom.r += v*dt;
    atom.v = v;
}

void Verlet::step_one_velocities(Atoms::iterator atomit, double dt) {
	auto& atom = *atomit;
	if (atom.tag == 0) return;
	vector a = atom.f; a /= atom.m;
    atom.v += a*0.5*dt;
}

} // namespace mdcraft::solver::stepper

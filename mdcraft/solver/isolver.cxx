#include <mdcraft/solver/isolver.h>


namespace mdcraft::solver {

void ISolver::evaluate(Atoms& atoms, NeibsList& nlist) {
    prepare(atoms, atoms, nlist);
    for (auto& stage : stages) stage(atoms, atoms, nlist);
}
#ifdef mdcraft_ENABLE_MPI
void ISolver::evaluate(Decomp& decomp, NeibsList& local, NeibsList& alien) {
    prepare(decomp.locals(), decomp.locals(), local);
    for (auto& stage : stages) {
        decomp.exchange_start();
        stage(decomp.locals(), decomp.locals(), local);
        decomp.exchange_end();
        stage(decomp.locals(), decomp.aliens(), alien);
    }
}
#endif

ISolver::~ISolver() = default;

void ISolver::prepare(Atoms& atoms, Atoms& neibs, NeibsList& nlist) { }

#ifdef mdcraft_ENABLE_MPI
void ISolver::prepare(Decomp& decomp, NeibsList& nlist1, NeibsList& nlist2) { }
#endif
} // namespace mdcraft::solver
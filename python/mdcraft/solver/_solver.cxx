#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
// #include <pybind11/stl.h>

#include <mdcraft/configuration.h>
#include <mdcraft/data/atom.h>
#include <mdcraft/tools/threads.h>
#include <mdcraft/neibs/verlet-list.h>
#include <mdcraft/solver/single.h>
#include <mdcraft/solver/multi.h>

#include <mdcraft/solver/dp-solver.h>

#ifdef mdcraft_ENABLE_MLIP
#include <mdcraft/solver/mlip.h>
#endif
#ifdef mdcraft_ENABLE_BPNN
#include <mdcraft/solver/ternary.h>
#endif
#ifdef mdcraft_ENABLE_MPI
#include <mdcraft/decomp/decomp.h>
#endif

namespace py = pybind11;

using mdcraft::tools::Threads;
using mdcraft::data::Atoms;

using NeibsList = mdcraft::neibs::List;
using mdcraft::solver::Single;
using mdcraft::solver::Multi;

#ifdef mdcraft_ENABLE_MLIP
using ::mdcraft::solver::MLIP;
using MTPRPotential = ::mdcraft::solver::potential::MTPR;
#endif

using Pair       = mdcraft::solver::Pair;
using Potential  = mdcraft::solver::potential::Base;
using Thermostat = mdcraft::solver::thermostat::Base;
using Boundary   = mdcraft::solver::boundary::Base;

using mdcraft::solver::ISolver;
using mdcraft::solver::ManyBodySolver;
using mdcraft::solver::potential::BaseManyBody;
#ifdef mdcraft_ENABLE_BPNN
using mdcraft::solver::Ternary;
using mdcraft::solver::potential::BaseTernary;
#endif
#ifdef mdcraft_ENABLE_MPI
using mdcraft::decomp::Decomp;
#endif

PYBIND11_MODULE(_mdcraft_solver, m) {

py::class_<ISolver>(m, "ISolver");

py::class_<Pair, ISolver>(m, "Pair")
	.def("prepare",   (void (Pair::*) (Atoms&, Atoms&, NeibsList&)) &Pair::prepare)
	.def("charge",    (void (Pair::*) (Atoms&, Atoms&, NeibsList&)) &Pair::charge)
	.def("forces",    (void (Pair::*) (Atoms&, Atoms&, NeibsList&)) &Pair::forces)
	.def("virials",   (void (Pair::*) (Atoms&, Atoms&, NeibsList&)) &Pair::virials)
#ifdef mdcraft_ENABLE_MPI
	.def("prepare",   (void (Pair::*) (Decomp&, NeibsList&, NeibsList&)) &Pair::prepare)
	.def("charge",    (void (Pair::*) (Decomp&, NeibsList&, NeibsList&)) &Pair::charge)
	.def("forces",    (void (Pair::*) (Decomp&, NeibsList&, NeibsList&)) &Pair::forces)
	.def("virials",   (void (Pair::*) (Decomp&, NeibsList&, NeibsList&)) &Pair::virials)
#endif
	;

py::class_<Single, Pair>(m, "Single")
	.def(py::init<
		Potential&,
		Boundary&,
		Thermostat&,
		Threads&
		>(), 
		py::arg("potential"),
		py::arg("boundary")   = mdcraft::solver::boundary::dummy_boundary,
		py::arg("thermostat") = mdcraft::solver::thermostat::dummy_thermostat,
		py::arg("threads")    = mdcraft::tools::dummy_pool,
        py::keep_alive<1, 2>(), py::keep_alive<1, 3>(),
        py::keep_alive<1, 4>(), py::keep_alive<1, 5>()
	);

py::class_<Multi, Pair>(m, "Multi")
	.def(py::init<Boundary&, Thermostat&, std::size_t, Threads&>(),
        py::arg("boundary") = mdcraft::solver::boundary::dummy_boundary,
        py::arg("thermostat") = mdcraft::solver::thermostat::dummy_thermostat,
        py::arg("nkinds") = 1,
        py::arg("threads") = mdcraft::tools::dummy_pool,
        py::keep_alive<1, 2>(), py::keep_alive<1, 3>(), py::keep_alive<1, 5>())
	.def("add_potential", [](Multi& multi, 
		std::size_t i,
		std::size_t j,
		Potential&  pot
	) {
		multi.add_potential(i,j,pot);
	}, py::keep_alive<1, 4>());

#ifdef mdcraft_ENABLE_MLIP
py::class_<MLIP, Pair>(m, "MLIP")
	.def(py::init<
		MTPRPotential&,
		Boundary&,
		Thermostat&,
		Threads&
	    >(),
		py::arg("potential"),
		py::arg("boundary")     = mdcraft::solver::boundary::dummy_boundary,
		py::arg("thermostat")   = mdcraft::solver::thermostat::dummy_thermostat,
		py::arg("threads")      = mdcraft::tools::dummy_pool
	)
	// .def("forces", &MLIP::forces,
	// 	py::arg("atoms"),
	// 	py::arg("neibs"),
	// 	py::arg("nlist"),
	// 	"Calculate forces for all atoms.")
	.def("virials", &MLIP::virials,
		py::arg("atoms"),
		py::arg("neibs"),
		py::arg("nlist"),
		"Calculate virials for all atoms."
	);
#endif


py::class_<ManyBodySolver, ISolver>(m, "ManyBodySolver")
	.def(py::init<
		BaseManyBody&,
		Thermostat&,
		Threads&
	>(),
		py::arg("potential"),
                py::arg("thermostat")   = mdcraft::solver::thermostat::dummy_thermostat,
                py::arg("threads")      = mdcraft::tools::dummy_pool
        )
	.def("prepare",   (void (ManyBodySolver::*) (Atoms&, Atoms&, NeibsList&)) &ManyBodySolver::prepare)
        .def("forces",  [](ManyBodySolver& self, Atoms& a, Atoms& aa, NeibsList& l) { return self.forces(a, aa, l); } );

#ifdef mdcraft_ENABLE_BPNN
py::class_<Ternary, ISolver>(m, "Ternary")
	.def(py::init<
		BaseTernary&,
		Thermostat&,
		Threads&
	>(),
		py::arg("potential"),
		py::arg("thermostat")   = mdcraft::solver::thermostat::dummy_thermostat,
		py::arg("threads")      = mdcraft::tools::dummy_pool
	)
	.def("prepare",   (void (Ternary::*) (Atoms&, Atoms&, NeibsList&)) &Ternary::prepare)
	.def("forces",  [](Ternary& self, Atoms& a, Atoms& aa, const NeibsList& l) { return self.forces(a, aa, l); } )
	.def("virials", [](Ternary& self, Atoms& a, Atoms& aa, const NeibsList& l) { return self.virials(a, aa, l); } );
#endif
}

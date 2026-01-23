#include <mdcraft/configuration.h>

#include <mdcraft/configuration.h>
#include <mdcraft/solver/potential/LJs.h>
#include <mdcraft/solver/potential/EAM.h>
#include <mdcraft/solver/potential/pair.h>
#include <mdcraft/tools/threads.h>

#ifdef mdcraft_ENABLE_CUDA
#include <mdcraft/solver/cuda/solver_cuda.h>
#endif

#ifdef mdcraft_ENABLE_BPNN
#include <mdcraft/solver/potential/maise_bpnn.h>
#endif

#ifdef mdcraft_ENABLE_DeePMD
#include <mdcraft/solver/potential/DeePMD.h>
#endif

#ifdef mdcraft_ENABLE_MLIP4
#include <mdcraft/solver/potential/mlip4.h>
#endif

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
// #include <pybind11/stl.h>

namespace py = pybind11;

using ::mdcraft::data::Atoms;
using NeibsList = ::mdcraft::neibs::List;

using ::mdcraft::tools::Threads;

using Potential = ::mdcraft::solver::potential::Base;
using LJs       = ::mdcraft::solver::potential::LJs;
using EAM       = ::mdcraft::solver::potential::EAM;
using Pair      = ::mdcraft::solver::potential::Pair;
#ifdef mdcraft_ENABLE_MLIP4
using MLIP4Potential        = ::mdcraft::solver::potential::MLIP4Potential;
#endif

using ::mdcraft::lattice::Domain;
using PotentialManyBody     = ::mdcraft::solver::potential::BaseManyBody;
#ifdef mdcraft_ENABLE_DeePMD
using DeepModelingPotential = ::mdcraft::solver::potential::DeepModelingPotential;
#endif

using PotentialTernary   = ::mdcraft::solver::potential::BaseTernary;
#ifdef mdcraft_ENABLE_BPNN
using MaiseBPNNPotential = ::mdcraft::solver::potential::MaiseBPNNPotential;
#endif

PYBIND11_MODULE(_mdcraft_solver_potential, m) {

py::class_<Potential>(m, "Potential")
	.def(py::init<double>(),
		 py::arg("Rcutoff") 
	)
	.def_property_readonly("rcut", [](Potential& pot) -> double {
		return pot.rcut();
	})
	;

py::class_<PotentialManyBody>(m, "PotentialManyBody")
	.def(py::init<double>(),
		 py::arg("Rcutoff") 
	)
	.def_property_readonly("rcut", [](PotentialManyBody& pot) -> double {
		return pot.rcut();
	})
	;

#ifdef mdcraft_ENABLE_BPNN
py::class_<PotentialTernary>(m, "PotentialTernary")
	.def(py::init<double>(),
		 py::arg("Rcutoff") 
	)
	.def_property_readonly("rcut", [](PotentialTernary& pot) -> double {
		return pot.rcut();
	})
	;
#endif

py::class_<LJs,Potential>(m, "LJs")
	.def(py::init<double, double, double>(),
		py::arg("aVr"),
		py::arg("rVr"),
		py::arg("Rcutoff")
	)
	.def("value", [](LJs& ljspot, py::array_t<double>& r2) -> py::array_t<double> {
		py::array_t<double> u(r2.size());
		auto cr2 = r2.data();
		auto cu = u.mutable_data();
		ljspot.value(cr2, cu, r2.size());
		return u;
	})
	.def("value", [](LJs& ljspot, double r) -> double {
		double u = 0.0;
		double r2 = r*r;
		ljspot.value(&r2, &u, 1);
		return u;
	});

py::class_<EAM,Potential>(m, "EAM")
	.def(py::init([](
			double Rcutoff,
            py::array_t<double> Vpar, 
            py::array_t<double> Rho, 
            py::array_t<double> Femb, 
            int    Nr1,
            double dr,
            int    Nrh1,
            double drh
        ) {
			// get c-array representation
            const double* cFemb = Femb.data();
            const double* cRho = Rho.data();
            const double* cVpar = Vpar.data();

            auto pot = new EAM(
                Rcutoff, cFemb, cRho, cVpar, Nr1, dr, Nrh1, drh
            );
            return pot;
        }),
			py::arg("Rcutoff"),
            py::arg("pairPotential"),
            py::arg("chargeDensity"),
            py::arg("embeddingEnergy"),
            py::arg("NrPot"),
            py::arg("drPot"),
            py::arg("NrEmb"),
            py::arg("drEmb")
        )
	.def(py::init([](std::string filename) {
			auto pot = new EAM(filename);
            return pot;
        }), 
		py::arg("filename")
	)
	// .def("U", [](EAM& potential, py::array_t<double>& r) -> py::array_t<double> {
	// 	py::array_t<double> u(r.size());
	// 	auto cr = r.data();
	// 	auto cu = u.mutable_data();
	// 	potential.U(cr, cu, r.size());
	// 	return u;
	// })
	.def("value", [](EAM& pot, double r) -> double {
		double u = 0.0;
		u = pot.value(::mdcraft::data::vector(r, 0, 0));
		return u;
	})
	.def("embeddingEnergy", [](EAM& pot, py::array_t<double>& r) -> py::array_t<double> {
		py::array_t<double> u(r.size());
		py::array_t<double> du(r.size());
		auto cr = r.data();
		auto cu = u.mutable_data();
		auto cdu = du.mutable_data();
		pot.embeddingEnergy(cr, cu, cdu, r.size());
		return u;
	});


py::class_<Pair,Potential>(m, "Pair")
	.def(py::init([](
			double Rcutoff,
            py::array_t<double> Vpar, 
            int    Nr1,
            double dr
        ) {
			// get c-array representation
            const double* cVpar = Vpar.data();

            auto pot = new Pair(
                Rcutoff, cVpar, Nr1, dr
            );
            return pot;
        }),
			py::arg("Rcutoff"),
            py::arg("pairPotential"),
            py::arg("NrPot"),
            py::arg("drPot")
        )
	.def(py::init([](std::string filename) {
			auto pot = new Pair(filename);
            return pot;
        }), 
		py::arg("filename")
	)
	.def("value", [](Pair& pot, double r) -> double {
		double u = 0.0;
		u = pot.value(::mdcraft::data::vector(r, 0, 0));
		return u;
	});

#ifdef mdcraft_ENABLE_DeePMD
	py::class_<DeepModelingPotential, PotentialManyBody>(m, "DeepModelingPotential")
		.def(py::init<const std::string&, Domain&>(),
		py::arg("pot_filename"),
        py::arg("domain")
    )
	.def("force", &DeepModelingPotential::force);
#endif

#ifdef mdcraft_ENABLE_MLIP4
	py::class_<MLIP4Potential, PotentialManyBody>(m, "MLIP4")
		.def(py::init<const std::string&, Domain&>(),
		py::arg("pot_filename"),
        py::arg("domain")
    )
	.def("force", &MLIP4Potential::force);
#endif

#ifdef mdcraft_ENABLE_BPNN
	py::class_<MaiseBPNNPotential, PotentialTernary>(m, "MaiseBPNNPotential")
		.def(py::init<std::string_view, Domain&>(),
		py::arg("pot_filename"),
        py::arg("domain")
    )
	.def("force", &MaiseBPNNPotential::force);
#endif

#ifdef mdcraft_ENABLE_CUDA_CRTP
py::class_<BasePotential>(m, "PotentialCUDA")
	.def(py::init([] (
			// ???
		) {
			return new BasePotential();
		})
	)
	;

py::class_<cuLJs, BasePotential, TBasePotential<cuLJs>::PtrHolder>(m, "cuLJs")
	.def(py::init([] (double a, double sigma, double rcut)
	{
		auto pot = TBasePotential<cuLJs>::create(a, sigma, rcut);
		return pot;
    }),
		py::arg("a"),
		py::arg("sigma"),
		py::arg("rcut")
	)
	.def_property("a", &cuLJs::get_a, &cuLJs::set_a)
	.def_property("rcut", &cuLJs::get_rcut, &cuLJs::set_rcut)
	.def_property("sigma", &cuLJs::get_sigma, &cuLJs::set_sigma)
	;

py::class_<cuEAM, BasePotential>(m, "cuEAM")
	.def(py::init<>())
	;
#endif

}

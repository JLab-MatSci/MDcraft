#include <mdcraft/configuration.h>
#include <mdcraft/solver/potential/LJs.h>
#include <mdcraft/solver/potential/LJs-cut.h>

#include <mdcraft/solver/potential/EAM.h>

#include <mdcraft/solver/potential/EAMW1046.h>
#include <mdcraft/solver/potential/EAMW0644.h>
#include <mdcraft/solver/potential/EAMV1622.h>

#include <mdcraft/solver/potential/EAMWMOD1046.h>
#include <mdcraft/solver/potential/EAMWMOD0644.h>
#include <mdcraft/solver/potential/EAMMOD1044.h>

#include <mdcraft/solver/potential/pair.h>
#include <mdcraft/tools/threads.h>

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
#include <cmath>
#include <pybind11/numpy.h>
// #include <pybind11/stl.h>

namespace py = pybind11;

using ::mdcraft::data::Atoms;
using NeibsList = ::mdcraft::neibs::List;

using ::mdcraft::tools::Threads;

using Potential = ::mdcraft::solver::potential::Base;
using LJs       = ::mdcraft::solver::potential::LJs;
using LJsCut       = ::mdcraft::solver::potential::LJsCut;

using EAM       = ::mdcraft::solver::potential::EAM;

using EAMW1046   = ::mdcraft::solver::potential::EAMW1046;
using EAMW0644   = ::mdcraft::solver::potential::EAMW0644;
using EAMV1622   = ::mdcraft::solver::potential::EAMV1622;

using EAMWMOD1046   = ::mdcraft::solver::potential::EAMWMOD1046;
using EAMWMOD0644   = ::mdcraft::solver::potential::EAMWMOD0644;
using EAMMOD1044 = ::mdcraft::solver::potential::EAMMOD1044;

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

static void require_array(const py::array_t<double>& a, py::ssize_t size) {
    if (a.ndim() != 1 || a.size() != size || a.strides(0) != sizeof(double))
        throw py::value_error("Expected equal-length contiguous 1D coefficient arrays");
    for (py::ssize_t i=0; i<a.size(); ++i)
        if (!std::isfinite(a.data()[i])) throw py::value_error("Coefficients must be finite");
}

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
        require_array(r2, r2.size());
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

py::class_<LJsCut,Potential>(m, "LJsCut")
	.def(py::init<double, double, double>(),
		py::arg("a"),
		py::arg("s"),
		py::arg("rcut")
	)
	;

m.attr("LJsTEST") = m.attr("LJsCut");

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
            require_array(Femb, Nrh1); require_array(Rho, Nr1); require_array(Vpar, Nr1);
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
	.def("value", [](EAM& pot, double r) -> double {
		double u = 0.0;
		u = pot.value(::mdcraft::data::vector(r, 0, 0));
		return u;
	})
	.def("embeddingEnergy", [](EAM& pot, py::array_t<double>& r) -> py::array_t<double> {
        require_array(r, r.size());
		py::array_t<double> u(r.size());
		py::array_t<double> du(r.size());
		auto cr = r.data();
		auto cu = u.mutable_data();
		auto cdu = du.mutable_data();
		pot.embeddingEnergy(cr, cu, cdu, r.size());
		return u;
	});


py::class_<EAMW1046, Potential>(m, "EAMW1046")
	.def(py::init([](
			double Rcutoff,
			double aVr,
			double bVr,
			double rVr,
            py::array_t<double> V,
            py::array_t<double> W,
            py::array_t<double> C,
            py::array_t<double> D,
            py::array_t<double> F,
            py::array_t<double> G
        ) {
			// get c-array representation
            const double* cV = V.data();
            const double* cW = W.data();
            const double* cC = C.data();
            const double* cD = D.data();
            const double* cF = F.data();
            const double* cG = G.data();

            std::size_t arrsz = V.size();
            require_array(V, arrsz); require_array(W, arrsz);
            require_array(C, arrsz); require_array(D, arrsz);
            require_array(F, arrsz); require_array(G, arrsz);

            auto pot = new EAMW1046(
                Rcutoff, aVr, bVr, rVr, arrsz, cV, cW, cC, cD, cF, cG
            );
            return pot;
        }),
			py::arg("Rcutoff"),
            py::arg("aVr"),
            py::arg("bVr"),
            py::arg("rVr"),
            py::arg("V"),
            py::arg("W"),
            py::arg("C"),
            py::arg("D"),
            py::arg("F"),
            py::arg("G")
        )

	.def("value", [](EAMW1046& pot, double r) -> double {
		double u = 0.0;
		u = pot.value(r);
		return u;
	})
	.def("value", [](EAMW1046& pot, py::array_t<double>& r2) -> py::array_t<double> {
        require_array(r2, r2.size());
		py::array_t<double> u(r2.size());
		auto cr2 = r2.data();
		auto cu = u.mutable_data();
		pot.value(cr2, cu, r2.size());
		return u;
	})
	.def("embeddingEnergy", [](EAMW1046& pot, py::array_t<double>& r) -> py::array_t<double> {
        require_array(r, r.size());
		py::array_t<double> u(r.size());
		py::array_t<double> du(r.size());
		auto cr = r.data();
		auto cu = u.mutable_data();
		auto cdu = du.mutable_data();
		pot.embeddingEnergy(cr, cu, cdu, r.size());
		return u;
	});


py::class_<EAMW0644, Potential>(m, "EAMW0644")
	.def(py::init([](
			double Rcutoff,
			double aVr,
			double bVr,
			double rVr,
            py::array_t<double> V,
            py::array_t<double> W,
            py::array_t<double> C,
            py::array_t<double> D,
            py::array_t<double> F,
            py::array_t<double> G
        ) {
			// get c-array representation
            const double* cV = V.data();
            const double* cW = W.data();
            const double* cC = C.data();
            const double* cD = D.data();
            const double* cF = F.data();
            const double* cG = G.data();

            std::size_t arrsz = V.size();
            require_array(V, arrsz); require_array(W, arrsz);
            require_array(C, arrsz); require_array(D, arrsz);
            require_array(F, arrsz); require_array(G, arrsz);

            auto pot = new EAMW0644(
                Rcutoff, aVr, bVr, rVr, arrsz, cV, cW, cC, cD, cF, cG
            );
            return pot;
        }),
			py::arg("Rcutoff"),
            py::arg("aVr"),
            py::arg("bVr"),
            py::arg("rVr"),
            py::arg("V"),
            py::arg("W"),
            py::arg("C"),
            py::arg("D"),
            py::arg("F"),
            py::arg("G")
        )

	.def("value", [](EAMW0644& pot, double r) -> double {
		double u = 0.0;
		u = pot.value(r);
		return u;
	})
	.def("value", [](EAMW0644& pot, py::array_t<double>& r2) -> py::array_t<double> {
        require_array(r2, r2.size());
		py::array_t<double> u(r2.size());
		auto cr2 = r2.data();
		auto cu = u.mutable_data();
		pot.value(cr2, cu, r2.size());
		return u;
	})
	.def("embeddingEnergy", [](EAMW0644& pot, py::array_t<double>& r) -> py::array_t<double> {
        require_array(r, r.size());
		py::array_t<double> u(r.size());
		py::array_t<double> du(r.size());
		auto cr = r.data();
		auto cu = u.mutable_data();
		auto cdu = du.mutable_data();
		pot.embeddingEnergy(cr, cu, cdu, r.size());
		return u;
	});

py::class_<EAMV1622, Potential>(m, "EAMV1622")
	.def(py::init([](
			double Rcutoff,
			double aVr,
			double bVr,
			double rVr,
            py::array_t<double> V,
            py::array_t<double> W
        ) {
			// get c-array representation
            const double* cV = V.data();
            const double* cW = W.data();

            std::size_t arrsz = V.size();
            require_array(V, arrsz); require_array(W, arrsz);


            auto pot = new EAMV1622(
                Rcutoff, aVr, bVr, rVr, arrsz, cV, cW
            );
            return pot;
        }),
			py::arg("Rcutoff"),
            py::arg("aVr"),
            py::arg("bVr"),
            py::arg("rVr"),
            py::arg("V"),
            py::arg("W")
        )

	.def("value", [](EAMV1622& pot, double r) -> double {
		double u = 0.0;
		u = pot.value(r);
		return u;
	})
	.def("value", [](EAMV1622& pot, py::array_t<double>& r2) -> py::array_t<double> {
        require_array(r2, r2.size());
		py::array_t<double> u(r2.size());
		auto cr2 = r2.data();
		auto cu = u.mutable_data();
		pot.value(cr2, cu, r2.size());
		return u;
	});

py::class_<EAMWMOD1046, Potential>(m, "EAMWMOD1046")
	.def(py::init([](
			double Rcutoff,
			double aVr,
			double bVr,
			double rVr,
            py::array_t<double> V,
            py::array_t<double> W,
            py::array_t<double> C,
            py::array_t<double> D,
            py::array_t<double> F,
            py::array_t<double> G
        ) {
			// get c-array representation
            const double* cV = V.data();
            const double* cW = W.data();
            const double* cC = C.data();
            const double* cD = D.data();
            const double* cF = F.data();
            const double* cG = G.data();

            std::size_t arrsz = V.size();
            require_array(V, arrsz); require_array(W, arrsz);
            require_array(C, arrsz); require_array(D, arrsz);
            require_array(F, arrsz); require_array(G, arrsz);

            auto pot = new EAMWMOD1046(
                Rcutoff, aVr, bVr, rVr, arrsz, cV, cW, cC, cD, cF, cG
            );
            return pot;
        }),
			py::arg("Rcutoff"),
            py::arg("aVr"),
            py::arg("bVr"),
            py::arg("rVr"),
            py::arg("V"),
            py::arg("W"),
            py::arg("C"),
            py::arg("D"),
            py::arg("F"),
            py::arg("G")
        )

	.def("value", [](EAMWMOD1046& pot, double r) -> double {
		double u = 0.0;
		u = pot.value(r);
		return u;
	})
	.def("value", [](EAMWMOD1046& pot, py::array_t<double>& r2) -> py::array_t<double> {
        require_array(r2, r2.size());
		py::array_t<double> u(r2.size());
		auto cr2 = r2.data();
		auto cu = u.mutable_data();
		pot.value(cr2, cu, r2.size());
		return u;
	})
	.def("embeddingEnergy", [](EAMWMOD1046& pot, py::array_t<double>& r) -> py::array_t<double> {
        require_array(r, r.size());
		py::array_t<double> u(r.size());
		py::array_t<double> du(r.size());
		auto cr = r.data();
		auto cu = u.mutable_data();
		auto cdu = du.mutable_data();
		pot.embeddingEnergy(cr, cu, cdu, r.size());
		return u;
	});


py::class_<EAMWMOD0644, Potential>(m, "EAMWMOD0644")
	.def(py::init([](
			double Rcutoff,
			double aVr,
			double bVr,
			double rVr,
            py::array_t<double> V,
            py::array_t<double> W,
            py::array_t<double> C,
            py::array_t<double> D,
            py::array_t<double> F,
            py::array_t<double> G
        ) {
			// get c-array representation
            const double* cV = V.data();
            const double* cW = W.data();
            const double* cC = C.data();
            const double* cD = D.data();
            const double* cF = F.data();
            const double* cG = G.data();

            std::size_t arrsz = V.size();
            require_array(V, arrsz); require_array(W, arrsz);
            require_array(C, arrsz); require_array(D, arrsz);
            require_array(F, arrsz); require_array(G, arrsz);

            auto pot = new EAMWMOD0644(
                Rcutoff, aVr, bVr, rVr, arrsz, cV, cW, cC, cD, cF, cG
            );
            return pot;
        }),
			py::arg("Rcutoff"),
            py::arg("aVr"),
            py::arg("bVr"),
            py::arg("rVr"),
            py::arg("V"),
            py::arg("W"),
            py::arg("C"),
            py::arg("D"),
            py::arg("F"),
            py::arg("G")
        )

	.def("value", [](EAMWMOD0644& pot, double r) -> double {
		double u = 0.0;
		u = pot.value(r);
		return u;
	})
	.def("value", [](EAMWMOD0644& pot, py::array_t<double>& r2) -> py::array_t<double> {
        require_array(r2, r2.size());
		py::array_t<double> u(r2.size());
		auto cr2 = r2.data();
		auto cu = u.mutable_data();
		pot.value(cr2, cu, r2.size());
		return u;
	})
	.def("embeddingEnergy", [](EAMWMOD0644& pot, py::array_t<double>& r) -> py::array_t<double> {
        require_array(r, r.size());
		py::array_t<double> u(r.size());
		py::array_t<double> du(r.size());
		auto cr = r.data();
		auto cu = u.mutable_data();
		auto cdu = du.mutable_data();
		pot.embeddingEnergy(cr, cu, cdu, r.size());
		return u;
	});

py::class_<EAMMOD1044, Potential>(m, "EAMMOD1044")
	.def(py::init([](
			double Rcutoff,
			double aVr,
			double bVr,
			double rVr,
            py::array_t<double> V,
            py::array_t<double> W,
            py::array_t<double> C,
            py::array_t<double> D,
            py::array_t<double> F,
            py::array_t<double> G
        ) {
			// get c-array representation
            const double* cV = V.data();
            const double* cW = W.data();
            const double* cC = C.data();
            const double* cD = D.data();
            const double* cF = F.data();
            const double* cG = G.data();

            std::size_t arrsz = V.size();
            require_array(V, arrsz); require_array(W, arrsz);
            require_array(C, arrsz); require_array(D, arrsz);
            require_array(F, arrsz); require_array(G, arrsz);

            auto pot = new EAMMOD1044(
                Rcutoff, aVr, bVr, rVr, arrsz, cV, cW, cC, cD, cF, cG
            );
            return pot;
        }),
			py::arg("Rcutoff"),
            py::arg("aVr"),
            py::arg("bVr"),
            py::arg("rVr"),
            py::arg("V"),
            py::arg("W"),
            py::arg("C"),
            py::arg("D"),
            py::arg("F"),
            py::arg("G")
        )
    .def("value", [](EAMMOD1044& pot, double r) { return pot.value(r); })
    .def("embeddingEnergy", [](EAMMOD1044& pot, double rho) {
        double energy; pot.embeddingEnergy(rho, energy); return energy;
    });


py::class_<Pair,Potential>(m, "Pair")
	.def(py::init([](
			double Rcutoff,
            py::array_t<double> Vpar, 
            int    Nr1,
            double dr
        ) {
			// get c-array representation
            require_array(Vpar, Nr1);
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

}

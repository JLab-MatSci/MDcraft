#include <iostream>
#include <utility>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <mpi4py/mpi4py.h>

#include <mdcraft/decomp/VD3.h>

namespace py = pybind11;

using mdcraft::lattice::Domain;
using mdcraft::data::Atoms;
using mdcraft::tools::Threads;

using mdcraft::decomp::Decomp;
using mdcraft::decomp::VD3;

std::vector<mdcraft::data::vector> convert_centers(const py::object& obj) {
	using mdcraft::data::vector;

	// 1. None object -> empty centers
	if (obj.is_none()) { return {}; }

	// 2. Если это list
	if (py::isinstance<py::list>(obj)) {
		auto lst = obj.cast<py::list>();
		py::size_t n = lst.size();
		std::vector<vector> result(n);
		for (py::size_t i = 0; i < n; ++i) {
			py::object item = lst[i];
			// 2a. Элемент — тоже list/tuple (например [x,y,z])
			if (py::isinstance<py::list>(item) || py::isinstance<py::tuple>(item)) {
				auto it = py::cast<py::list>(item);
				if (it.size() != 3) {
					throw std::runtime_error("Each center must have 3 coords");
				}
				double x = py::cast<double>(it[0]);
				double y = py::cast<double>(it[1]);
				double z = py::cast<double>(it[2]);
				result[i] = {x, y, z};
				continue;
			}
			// 2b. Элемент — numpy.array([x,y,z])
			if (py::isinstance<py::array>(item)) {
				auto arr = item.cast<py::array>();
				if (arr.ndim() != 1 || arr.shape(0) != 3) {
					throw std::runtime_error("Each center must have 3 coords");
				}
				auto buf = arr.request();
				auto ptr = static_cast<double*>(buf.ptr);
				result[i] = {ptr[0], ptr[1], ptr[2]};
				continue;
			}
			throw std::runtime_error("Unsupported center type");
		}
		return result;
	}

    // 2. Если это numpy.ndarray
	if (py::isinstance<py::array>(obj)) {
        auto arr = obj.cast<py::array>();
        if (arr.ndim() != 2 || arr.shape(1) != 3) {
            throw std::runtime_error("Expected array of shape (N, 3)");
        }
        auto buf = arr.request();
        auto ptr = static_cast<double*>(buf.ptr);
        auto n = arr.shape(0);

    	std::vector<vector> result(n);
        for (py::size_t i = 0; i < n; ++i) {
            result[i] = {ptr[0], ptr[1], ptr[2]};
            ptr += 3;
        }
        return result;
    }

    throw std::runtime_error("Unsupported centers type, expected list of triplets or (N, 3) array)");
}

PYBIND11_MODULE(_mdcraft_decomp, m) {

py::class_<Decomp>(m, "Decomp")
	.def_property_readonly("rank", [](Decomp& D) -> int {
		return D.rank();
	})
    .def_property_readonly("size", [](Decomp& D) -> int {
        return D.size();
	})
	.def_property_readonly("locals", [](Decomp& D) -> Atoms& {
		return D.locals();
	}, py::return_value_policy::reference_internal)
	.def_property_readonly("aliens", [](Decomp& D) -> Atoms& {
		return D.aliens();
	}, py::return_value_policy::reference_internal)
	.def_property_readonly("border", [](Decomp& D) -> Atoms& {
		return D.border();
	}, py::return_value_policy::reference_internal)
	.def_property("measurer", nullptr, [](Decomp& D, const std::string& type) {
		D.set_measurer(type);
	})
	.def("exchange", [](Decomp& D) { D.exchange(); })
	.def("exchange_start", [](Decomp& D) { D.exchange_start(); })
	.def("exchange_end",   [](Decomp& D) { D.exchange_end(); })
	.def("update", [](Decomp& D, bool verbose) {
		D.update(verbose);
	},
		py::arg("verbose") = false
	)
	.def("prebalancing",   [](Decomp& D, int n_iters, bool verbose) {
		D.prebalancing(n_iters, verbose);
	},
		py::arg("n_iters") = 15,
		py::arg("verbose") = false
	)
;



py::class_<VD3, Decomp>(m, "VD3")
    .def(py::init([](
            py::object         comm_obj,
            Atoms&             atoms,
            Domain&            domain,
            int                dimension,
            py::object         centers,
            Threads&           pool,
            double             mobility,
            double             centroidal,
            double             growth_rate,
            const std::string& measurer
    ) {
    	// extract communicator from comm_obj
    	MPI_Comm comm = ((PyMPICommObject*) comm_obj.ptr())->ob_mpi;
    	// fill coords if needed
	    auto cpp_centers = convert_centers(centers);

    	VD3* vd3 = new VD3(
       		comm,
		    atoms,
			domain,
			pool,
			{
				.dimension   = dimension,
        		.mobility    = mobility,
        		.centroidal  = centroidal,
        		.growth_rate = growth_rate
			},
			cpp_centers
        );
    	vd3->set_measurer(measurer);
    	return vd3;
    }), py::arg("comm"), 
        py::arg("atoms"), 
        py::arg("domain"),
        py::arg("dimension") = 1,
        py::arg("centers") = py::none(),
        py::arg("threads") = mdcraft::tools::dummy_pool,
        py::arg("mobility") = 0.2,
        py::arg("centroidal") = 0.25,
        py::arg("growth_rate") = 0.02,
        py::arg("measurer") = "time"
    );

}

#include <mdcraft/configuration.h>
#include <pybind11/pybind11.h>
#include <mdcraft/tools/threads.h>
#ifdef mdcraft_ENABLE_MPI
#include <mpi4py/mpi4py.h>
#endif

namespace py = pybind11;

using mdcraft::tools::Threads;

PYBIND11_MODULE(_mdcraft_tools, m) {

py::class_<Threads>(m, "Threads")
    .def(py::init([](int count, py::object comm_obj) {
#ifndef mdcraft_ENABLE_MPI
        auto pool = new Threads(count);
        return pool;
#else
        if (comm_obj.is_none()) {
            auto pool = new Threads(count);
            return pool;
        }
        else {
            MPI_Comm comm = ((PyMPICommObject*) comm_obj.ptr())->ob_mpi;
            auto pool = new Threads(count, comm);
            return pool;
        }
#endif
    }),
        py::arg("count") = 100000,    // all available threads
        py::arg("comm")  = py::none() // MPI_COMM_WORLD by default
    )
    .def_property_readonly("count", [](Threads& pool) -> int {
        return pool.count();
    })
    .def("on",     [](Threads& pool, int count) { pool.on(); })
    .def("off",    [](Threads& pool) { pool.off(); });
}

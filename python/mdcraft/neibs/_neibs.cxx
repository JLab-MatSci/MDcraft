#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <mdcraft/data/atom.h>
#include <mdcraft/neibs/verlet-list.h>
#include <mdcraft/tools/threads.h>

namespace py = pybind11;

using mdcraft::lattice::Domain;

using mdcraft::data::Atoms;
using mdcraft::tools::Threads;
using NeibsList = mdcraft::neibs::List;
using mdcraft::neibs::VerletList;

PYBIND11_MODULE(_mdcraft_neibs, m) {

py::class_<NeibsList>(m, "NeibsList");

py::class_<VerletList>(m, "VerletList", "Verlet List class")
    .def(py::init([](
        Atoms&   atoms,
        Atoms&   neibs,
        Domain&  domain,
        Threads& threads,
        bool     half
    ) {
        auto nlist = new VerletList(
            atoms,
            neibs,
            domain,
            threads,
            half
        );
        return nlist;
    }),
        py::arg("atoms"),
        py::arg("neighbors"),
        py::arg("domain"),
        py::arg("threads")    = mdcraft::tools::dummy_pool,
        py::arg("half")       = false,
        R"(class constructor

        Parameters
        ----------

        atoms : Atoms, required
            Atoms for which list of neighbors will be calculated
        neighbors : Atoms, required
            Atoms among which neighbors will be searched
        dimension : int, optional
            Problem dimension
        threads : Threads, optional
            Threads. See User guide for more information.
        half : bool, optional
            Is this a half- or full-list?

    )")
    .def("update", [](VerletList& nlist,
        double kbuf
    ) -> void {
        nlist.update(kbuf);
    },
        py::arg("kbuf") = -1.0,
        R"(Updates list of neighbors

        Parameters
        ----------

        kbuf : double, optional
             Determines how many times the radius of the 
             neighbor search is greater than the doubled
             maximum atom size. If omitted, neibsSearchRadius
             is not updated.

    )")
    .def("sort", [](
        VerletList& nlist,
        Atoms&  elems
    ) -> Atoms {
        return nlist.sort(elems);
    },
        R"(Provides internal sorting for acceleration

        Parameters
        ----------

        elems : Atoms, required
            Atoms to be sorted.

        Returns
        ----------

        out : Atoms
            Sorted atoms

    )")
    .def("get", [](VerletList& nlist) -> NeibsList& {
        return nlist.get();
    }, py::return_value_policy::reference_internal)
    .def_property_readonly("half", [](VerletList& nlist) -> bool {
        return nlist.is_half();
    })
    .def("decide", [](VerletList& nlist) -> bool {
        return nlist.decide();
    },
        R"(Checks whether the current neighbor list should be rebuilt

        Returns
        -------

        out : bool
            True if at least one atom moved more than half of the
            current buffer distance since the last update().

    )")
    .def("__getitem__", [](VerletList& nlist, std::size_t i) -> py::array {
        auto result = nlist[i];
        auto ptr = result.data();
        std::vector<std::size_t> shape = {result.size()};
        std::vector<std::size_t> strides = {sizeof(::mdcraft::neibs::point_id)};
        return py::array_t<::mdcraft::neibs::point_id>(shape, strides, ptr);
    })
    .def("list_for", [](VerletList& nlist, double x) -> py::array {
        auto result = nlist.list_for({x, 0., 0.});
        auto np_result = py::array_t<std::size_t>(result.size());
        std::memcpy(np_result.mutable_data(), result.data(), sizeof(std::size_t)*result.size());
        return np_result;
    },
        R"(Returns list of neighbors for specified point

        Parameters
        ----------

        x : double, required
             x-coordinate

        Returns
        -------

        out : ndarray
            An array object with neighbors

    )")
    .def("list_for", [](VerletList& nlist, double x, double y) -> py::array {
        auto result = nlist.list_for({x, y, 0.});
        auto np_result = py::array_t<std::size_t>(result.size());
        std::memcpy(np_result.mutable_data(), result.data(), sizeof(std::size_t)*result.size());
        return np_result;
    },
        R"(Returns list of neighbors for specified atoms

        Parameters
        ----------

        x : double, required
             x-coordinate
        y : double, required
             y-coordinate

        Returns
        ----------

        out : ndarray
            An array object with neighbors

    )")
    .def("list_for", [](VerletList& nlist, double x, double y, double z) -> py::array {
        auto result = nlist.list_for({x, y, z});
        auto np_result = py::array_t<std::size_t>(result.size());
        std::memcpy(np_result.mutable_data(), result.data(), sizeof(std::size_t)*result.size());
        return np_result;
    },
        R"(Returns list of neighbors for specified atoms

        Parameters
        ----------

        x : double, required
             x-coordinate
        y : double, required
             y-coordinate
        z : double, required
             z-coordinate

        Returns
        ----------
        
        out : ndarray
            An array object with neighbors

    )")
    ;
}

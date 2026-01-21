#include <mdcraft/configuration.h>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <mdcraft/data/atom.h>

namespace py = ::pybind11;
using ::mdcraft::data::Atoms;

PYBIND11_MODULE(_mdcraft_data, m) {

struct Vector3D {
    double x;
    double y;
    double z;
};

struct Matrix3D {
    double xx;
    double xy;
    double xz;
    double yx;
    double yy;
    double yz;
    double zx;
    double zy;
    double zz;
};

struct Atom {
    Vector3D       r;
    Vector3D       v;
    Vector3D       f;
    Matrix3D       V;
    double         m;
    double         rcut;
    double         rns;
    double         Ep;
    double         Em;
    double         nc;
    unsigned short bc;
    unsigned short kind;
    unsigned int   uid;
    unsigned short tag = std::numeric_limits<unsigned short>::max();
};

PYBIND11_NUMPY_DTYPE(Vector3D,
    x,
    y,
    z
);

PYBIND11_NUMPY_DTYPE(Matrix3D,
    xx,
    xy,
    xz,
    yx,
    yy,
    yz,
    zx,
    zy,
    zz
);

PYBIND11_NUMPY_DTYPE(Atom, 
	r,
    v,
    f,
    V,
    m,
    rcut,
    rns,
    Ep,
    Em,
    nc,
    bc,
    kind,
    uid,
    tag
);

auto& api = py::detail::npy_api::get();
auto dtype = py::dtype::of<Atom>();

py::class_<Atoms>(m, "Atoms", py::buffer_protocol(), 
    "Atoms container class"
    )
    .def(py::init([](std::size_t size) {
        auto atoms = new Atoms(size);
        return atoms;
    }),
        py::arg("size") = 0,
        R"(class constructor

        Parameters
        ----------
        size : int, optional
            The number of atoms to be stored
        )"
    )
    .def_property_readonly("data", [api, dtype](Atoms& atoms) -> py::array {
        std::size_t atoms_shape[1] = {atoms.size()};
        std::size_t atoms_strides[1] = {sizeof(Atom)};
        dtype.inc_ref();
        return py::reinterpret_steal<py::array>(
            api.PyArray_NewFromDescr_(
                api.PyArray_Type_, 
                dtype.ptr(),
                1, 
                (Py_intptr_t*) atoms_shape, 
                (Py_intptr_t*) atoms_strides,
                (void *)       atoms.data(), 
                py::detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_ |
                py::detail::npy_api::NPY_ARRAY_ALIGNED_ |
                py::detail::npy_api::NPY_ARRAY_WRITEABLE_,
                nullptr
            )
        );
    },
        R"(Allows to access the underlying data as a NumPy array
            
            Returns
            -------

            NumPy array which points to the inner data of atoms container.
            It allows to modify the inner data by copying data from outer
            NumPy arrays. To list the avaliable fields of this structured array
            use command: print(atoms.data.dtype)
        )"
    )
     .def("resize", [](Atoms& atoms, std::size_t size) -> void {
        atoms.resize(size);
    },
        R"(Allows to change the size of underlying array

            Parameters
            ----------

            size : int
                New size of the underlying data storage.
            
        )"
    )
    .def("__len__", [](Atoms& atoms) -> int {
        return atoms.size();
    })
    .def("erase", [](Atoms& atoms, py::array_t<long long int> index) -> void {
        auto        index_data      = index.data();
        std::size_t index_count     = index.size();
        auto        atoms_data  = atoms.data();
        std::size_t atoms_count = atoms.size();
        for (std::size_t i = 0; i < index_count; i++) {
            std::size_t chunk_size = i + 1 == index_count ?
                                     atoms_count : index_data[i + 1];
            chunk_size -= index_data[i] + 1;
            if (index_data[i] + 1 < atoms_count)
                std::memmove(
                    (atoms_data + index_data[i] - i),
                    (atoms_data + index_data[i] + 1),
                    sizeof(Atom) * chunk_size
                );
        }
        atoms.resize(atoms_count - index_count);
    },
        R"(Allows to remove atoms according to the index list 

            Parameters
            ----------

            index : NumPy array
                The array of integer indices of atoms to remove.
            
        )"
    )
    .def("erase", [](Atoms& atoms, py::tuple index_tuple) -> void {
        py::array_t<long long int> index(index_tuple[0]);
        auto        index_data      = index.data();
        std::size_t index_count     = index.size();
        auto        atoms_data  = atoms.data();
        std::size_t atoms_count = atoms.size();
        for (std::size_t i = 0; i < index_count; i++) {
            std::size_t chunk_size = i + 1 == index_count ?
                                     atoms_count : index_data[i + 1];
            chunk_size -= index_data[i] + 1;
            if (index_data[i] + 1 < atoms_count)
                std::memmove(
                    (atoms_data + index_data[i] - i),
                    (atoms_data + index_data[i] + 1),
                    sizeof(Atom) * chunk_size
                );
        }
        atoms.resize(atoms_count - index_count);
    },
        R"(Allows to remove atoms according to the index list 

            Parameters
            ----------

            index : tuple
                The tuple, first element of which is the array of integer
                indices of atoms to remove. Usually returned by commands:
                - index = numpy.where(condition)
                - index = numpy.nonzero(elements)
        )"
    )
    .def("erase", [](Atoms& atoms, std::size_t from, std::size_t to) -> void {
        auto pfrom = atoms.begin() + from;
        auto pto   = atoms.begin() + to;
        atoms.erase(pfrom, pto);
    },
        R"(Allows to remove atoms range

            Parameters
            ----------

            from : int
                The index, from which the atoms are removed
            to : int
                The index, to which the atoms are removed
        )"
    )
   .def("join", [](Atoms& atoms, Atoms& other) -> void {
        atoms.insert(atoms.end(), other.begin(), other.end());
    },
        R"(Allows to join one atom storage to another

            Parameters
            ----------

            other : Atoms
                The atoms container, which is joined to the caller
        )"
    )
    ;

}
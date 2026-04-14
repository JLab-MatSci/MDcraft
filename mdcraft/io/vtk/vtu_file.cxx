#include <fstream>
#include <utility>

#include <mdcraft/io/vtk/vtu_file.h>
#include <mdcraft/tools/filesystem.h>

namespace mdcraft::io {

namespace fs = mdcraft::tools::filesystem;

using namespace data;

/// ===========================================================================
///             Реализация записи в виде набора статических функций
///             Implementation of writing as a set of static functions
/// ===========================================================================

inline bool is_big_endian() {
    union {
        uint32_t i;
        char c[4];
    } bint = {0x01020304};

    return bint.c[0] == 1;
}

inline std::string byteorder() {
    return is_big_endian() ? "BigEndian" : "LittleEndian";
}

namespace {
using datasize_t = std::uint32_t;  // The size of data buffer (do not change due to VTK)
using offset_t   = std::uint64_t;  // The type of arrays shifts in bytes
using index_t    = std::uint64_t;  // The type of primitives' indices
using type_t     = std::uint8_t;   // The primitive type

using byte_ptr = char*;

void write_header(std::ofstream &file, const Atoms &atoms, const Variables &variables) {
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" + byteorder() + "\">\n";
    file << "  <UnstructuredGrid>" << '\n';
    file << "    <Piece NumberOfPoints=\"" << atoms.size() << "\" NumberOfCells=\"" << atoms.size() << "\">\n";

    // Points
    offset_t byte_offset = 0;
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
    file << "      </Points>\n";
    byte_offset += sizeof(datasize_t) + 3 * atoms.size() * sizeof(double);

    // Cells
    file << "      <Cells>" << '\n';
    file << "        <DataArray type=\"" << VtkType::get<index_t>() << "\" Name=\"connectivity\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
    byte_offset += sizeof(datasize_t) + atoms.size() * sizeof(index_t);

    file << "        <DataArray type=\"" << VtkType::get<index_t>() << "\" Name=\"offsets\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
    byte_offset += sizeof(datasize_t) + atoms.size() * sizeof(index_t);

    file << "        <DataArray type=\"" << VtkType::get<type_t>() << "\" Name=\"types\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
    byte_offset += sizeof(datasize_t) + atoms.size() * sizeof(type_t);
    file << "      </Cells>\n";

    // Points Data
    file << "      <PointData>\n";
    for (auto &field: variables.list()) {
        file << "        <DataArray type=\"" << field.type() << "\" Name=\"" << field.name();
        if (!field.is_scalar()) {
            file << "\" NumberOfComponents=\"" << field.n_components();
        }
        file << "\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
        byte_offset += sizeof(datasize_t) + atoms.size() * field.size();
    }
    file << "      </PointData>\n";

    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
}

void write_primitives(std::ofstream &file, const Atoms &atoms) {
    // AppendedData
    file << "  <AppendedData encoding=\"raw\">\n";
    file << "_";

    // Points
    datasize_t data_size = 3 * atoms.size() * sizeof(double);
    file.write((byte_ptr) &data_size, sizeof(datasize_t));

    for (auto &atom: atoms) {
        file.write((byte_ptr) atom.r.data(), 3 * sizeof(double));
    }

    // Cells
    // Connectivity
    data_size = atoms.size() * sizeof(index_t);
    file.write((byte_ptr) &data_size, sizeof(datasize_t));

    for (index_t i = 0; i < atoms.size(); ++i) {
        file.write((byte_ptr) &i, sizeof(index_t));
    }

    // Offsets
    data_size = atoms.size() * sizeof(index_t);
    file.write((byte_ptr) &data_size, sizeof(datasize_t));

    for (size_t i = 1; i <= atoms.size(); ++i) {
        file.write((byte_ptr) &i, sizeof(index_t));
    }

    // Types
    data_size = atoms.size() * sizeof(type_t);
    file.write((byte_ptr) &data_size, sizeof(datasize_t));

    uint8_t type = 1;
    for (index_t i = 0; i < atoms.size(); ++i) {
        file.write((byte_ptr) &type, sizeof(type_t));
    }
}

void write_data(std::ofstream &file, Atoms &atoms, const Variables &variables) {
    std::vector<std::byte> temp;

    for (auto &field: variables.list()) {
        size_t field_size = field.size();
        datasize_t data_size = atoms.size() * field_size;

        file.write((byte_ptr) &data_size, sizeof(datasize_t));

        temp.resize(data_size);

        size_t counter = 0;
        for (auto& atom: atoms) {
            field.write(atom, temp.data() + counter * field_size);
            ++counter;
        }

        file.write((byte_ptr) temp.data(), data_size);
    }
}
}

/// ===========================================================================
///                            Функции класса
///                           The class methods
/// ===========================================================================

VtuFile::VtuFile(std::string filename) :
    filename(std::move(filename)) {
}

void VtuFile::save(Atoms &atoms, const Variables& variables) const {
    save(filename, atoms, variables);
}

void VtuFile::save(std::string filename, Atoms &atoms, const Variables &variables) {
    // Add file extension if it's absent
    if (filename.size() <= 4) {
        filename += ".vtu";
    }
    else {
        if (filename.substr(filename.size() - 4) != ".vtu") {
            filename += ".vtu";
        }
    }

    // Create folders if a complex filename given
    fs::path file_path(filename);
    fs::path dir_path = file_path.parent_path();

    if (!dir_path.empty()) {
        fs::create_directories(dir_path);
    }

    std::ofstream file(filename, std::ios::out | std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Warning: Cannot open file '" << filename << "'\n";
        return;
    }

    write_header(file, atoms, variables);
    write_primitives(file, atoms);
    write_data(file, atoms, variables);

    file.close();
}

} // namespace mdcraft::io

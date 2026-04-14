#include <iomanip>
#include <fstream>

#include <mdcraft/tools/network.h>
#include <mdcraft/tools/filesystem.h>
#include <mdcraft/io/vtk/pvd_file.h>
#include <mdcraft/io/vtk/vtu_file.h>

namespace mdcraft::io {

using namespace data;
using tools::Network;

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

PvdFile::PvdFile() : m_open(false), m_counter(0) { }

PvdFile::PvdFile(const std::string& filename, const std::string& directory)
    : PvdFile() {
#ifdef mdcraft_ENABLE_MPI
    m_net = Network(MPI_COMM_WORLD);
    open(filename, directory, !m_net.single());
#else
    open(filename, directory, false);
#endif
}

#ifdef mdcraft_ENABLE_MPI
PvdFile::PvdFile(Network network, const std::string& filename, const std::string& directory)
    : PvdFile() {
    m_net = network;
    open(filename, directory, !m_net.single());
}
#endif

void PvdFile::open(const std::string& filename, const std::string& _directory, bool distributed) {
    namespace fs = mdcraft::tools::filesystem;

    if (m_open) { return; }

    fs::path directory = fs::current_path();
    if (!_directory.empty()) {
        fs::path dir = _directory;
        if (dir.is_relative()) {
            directory /= dir;
        } else {
            directory = dir;
        }
    }

    // Мастер проверяет наличие директории и создает при необходимости
    // Master process checks if the directory exists and creates if needed
    if (m_net.master()) {
        if (!fs::exists(directory) || !fs::is_directory(directory)) {
            fs::create_directories(directory);
        }
    }

    m_filename = filename;
    if (filename.size() > 4) {
        if (filename.substr(filename.size() - 4) == ".pvd") {
            m_filename = filename.substr(filename.size() - 4);
        }
    }
    m_fullname = (directory / filename).string();

    // Мастер-процесс пишет заголовок PVD
    // Master-process writes a PVD-header
    m_distributed = m_net.single() ? false : distributed;
    if (m_distributed && !m_net.master()) {
        return;
    }

    /// Откроем файл и запишем заголовок
    // Open the file and write the header
    std::ofstream ofs;
    ofs.open(m_fullname + ".pvd");

    if (!ofs.is_open()) {
        std::cerr << "Warning: Cannot open .pvd file " << m_fullname << ".pvd\n";
        return;
    }

    ofs << "<?xml version=\"1.0\"?>\n";
    ofs << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"" + byteorder() + "\">\n";
    ofs << "    <Collection>" << std::endl;

    m_pos = ofs.tellp();

    ofs << "    </Collection>\n";
    ofs << "</VTKFile>\n";

    ofs.close();

    m_open = true;
}

void PvdFile::save(Atoms& atoms, double timestep) {
    VtuFile::save(get_filename(), atoms, variables);
    update_pvd(timestep);
}

std::string PvdFile::get_filename() const {
    std::string filename = m_fullname + "_" + std::to_string(m_counter);

    if (m_distributed) {
        filename += ".pt" + std::to_string(m_net.rank());
    }

    filename += ".vtu";
    return filename;
}

void PvdFile::update_pvd(double timestep) {
    // Мастер-процесс пишет PVD
    // Master-process writes to a PVD-file
    if (m_distributed && !m_net.master()) {
        ++m_counter;
        return;
    }

    if (!m_open) {
        std::string message = "PvdFile::save() error: You need to open PvdFile";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }

    std::fstream ofs;
    ofs.open(m_fullname + ".pvd");

    if (!ofs.is_open()) {
        std::cerr << "Cannot open file " << m_fullname << ".pvd\n";
    }

    ofs.seekg(m_pos, std::ios::beg);

    ofs << std::scientific << std::setprecision(15);

    if (m_distributed) {
        for (int r = 0; r < m_net.size(); ++r) {
            ofs << "        <DataSet timestep=\"" << timestep << "\" part=\"" << r << "\" file=\""
                << m_filename << "_" << m_counter << ".pt" << r << ".vtu" << "\"/>\n";
        }
    }
    else {
        ofs << "        <DataSet timestep=\"" << timestep << "\" part=\"0\" file=\""
            << m_filename << "_" << m_counter << ".vtu" << "\"/>\n";
    }

    m_pos = ofs.tellg();

    ofs << "    </Collection>\n";
    ofs << "</VTKFile>\n";

    ofs.close();

    ++m_counter;
}

} // namespace mdcraft::io

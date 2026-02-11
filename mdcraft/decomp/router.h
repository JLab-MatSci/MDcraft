#pragma once
#include <vector>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <future>
#include <mpi.h>
#include <cassert>
#include <mdcraft/data/atom.h>
#include <mdcraft/tools/network.h>

namespace mdcraft::decomp {

/// @brief Array of MPI isend/irecv requests
class Requests {
public:
    Requests() = default;
    /// @brief Array of 'null' requests
    explicit Requests(int size)
        : m_requests(size, MPI_REQUEST_NULL) { }
    /// @brief Get request by rank
    MPI_Request& operator[](int r) { return m_requests[r]; }
    /// @brief Wait for completion of all requests
    void wait();
private:
    std::vector<MPI_Request> m_requests;
};

/// @brief Manages communication operations
class Router {
public:
    /// @brief Arrays are initialized with zeros
    explicit Router(MPI_Comm comm);
    /// @brief Deletes MPI_Atom data type
    ~Router();
    /// @brief Set the number of elements to send
    void set_send_count(const std::vector<size_t>& send_count);
    /// @brief Set the number of elements to send
    void set_send_count(const std::atomic_size_t* send_count, int size);
    /// @brief Set the number of elements to receive
    void set_recv_count(const std::vector<size_t>& recv_count);
    /// @brief Assemble transfers
    void fill_partial();
    /// @brief Assemble complete transfer matrix
    void fill_complete();
    /// @brief Number of processes
    int size() const { return m_comm.size(); }
    /// @brief Is complete transfer matrix available?
    bool complete() const { return !m_send_recv.empty(); }
    /// @brief Number of transfers from process i to process j
    size_t get(int i, int j) const;
    /// @brief Number of transfers from process i to process j
    size_t operator()(int i, int j) const { return get(i, j); }
    /// @brief Required buffer size for sending messages
    size_t send_buffer_size() const;
    /// @brief Required buffer size for receiving messages
    size_t recv_buffer_size() const;
    const std::vector<size_t>& send_count() const { return m_send_count; }
    const std::vector<size_t>& recv_count() const { return m_recv_count; }
    const std::vector<size_t>& send_offset() const { return m_send_offset; }
    const std::vector<size_t>& recv_offset() const { return m_recv_offset; }
    size_t send_count(int r) const { return m_send_count[r]; }
    size_t recv_count(int r) const { return m_recv_count[r]; }
    /// @brief Print transfer information
    void print() const;
    /// @brief Asynchronous send
    Requests isend(const data::Atoms& src);
    /// @brief Asynchronous receive
    Requests irecv(data::Atoms& dst);
protected:
    /// @brief Print partial transfer matrix to console
    void print_partial() const;
    /// @brief Print complete transfer matrix to console
    void print_complete() const;
    tools::Network m_comm;  ///< MPI communicator
    MPI_Datatype MPI_Atom;  ///< MPI data type for Atom
    std::vector<size_t> m_send_count;   ///< Number of elements to send
    std::vector<size_t> m_recv_count;   ///< Number of elements to receive
    std::vector<size_t> m_send_offset;  ///< Offsets in send array
    std::vector<size_t> m_recv_offset;  ///< Offsets in receive array
    /// @brief Complete transfer matrix (optional)
    std::vector<size_t> m_send_recv;
};
} // namespace mdcraft::decomp

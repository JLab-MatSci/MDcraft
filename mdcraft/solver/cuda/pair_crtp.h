#pragma once

#include <numeric>

#include "isolver.h"

namespace mdcraft::solver::cuda {

template <typename P>
using Potential = mdcraft::solver::potential::cuda::TBasePotential<P>;

class PairSolver : public TBaseSolver<PairSolver>
{
public:
	PairSolver() = default;
};

template <typename S, typename P>
class TPairSolver : public PairSolver
{
public:
	using solver_t    = S;
	using potential_t = P;

	TPairSolver() = default;

#ifdef mdcraft_ENABLE_CUDA_THRUST
	void forces(Atoms& atoms, Atoms& neibs, NeibsList& nlist) /*const*/
	{
		static_assert(std::is_base_of_v<TPairSolver<solver_t, potential_t>, solver_t>);
		static_assert(std::is_trivially_copyable_v<std::remove_reference_t<solver_t>>);
		static_assert(std::is_standard_layout_v<std::remove_reference_t<solver_t>>);

		auto n_atoms = atoms.size();

		// convert to AtomsDev structure
		AtomsDev atoms_dev(atoms.begin(), atoms.end());

		auto first = thrust::make_zip_iterator(
			thrust::make_tuple(thrust::counting_iterator<std::size_t>(0), atoms_dev.begin())
		);
		auto last  = thrust::make_zip_iterator(
			thrust::make_tuple(thrust::counting_iterator<std::size_t>(n_atoms), atoms_dev.end())
		);

		DISPATCH_SOLVER_FOR_CUDA(solver_t, p_solver_dev);

		// get raw data under device Atoms verion
		data::Atom* p_atoms_dev = thrust::raw_pointer_cast(atoms_dev.data());

		// evaluate summed-up size of all the NlistOnes
		auto flatten_size = std::transform_reduce(nlist.begin(), nlist.end(), 0ull,
			std::plus<>(), [] (auto&& l) { return l.size(); }
		);
		auto&& [v_flat_dev, v_lens_dev, v_offs_dev] = neibs::create_nlist_device(
			nlist.begin(), nlist.end(), flatten_size
		);
		auto* p_flat_dev = thrust::raw_pointer_cast(v_flat_dev.data());
		auto* p_lens_dev = thrust::raw_pointer_cast(v_lens_dev.data());
		auto* p_offs_dev = thrust::raw_pointer_cast(v_offs_dev.data());

		thrust::for_each(first, last,
			[p_solver_dev, p_atoms_dev,
			 p_flat_dev, p_lens_dev, p_offs_dev] __device__ (auto&& t)
			{
				auto i = thrust::get<0>(t);
				auto& atom = thrust::get<1>(t);

				neibs::NeibsListOneDev nlist_dev
				{
					p_flat_dev + p_offs_dev[i],
					p_lens_dev[i]
				};

				// evaluate on device
				p_solver_dev->force_one(atom, p_atoms_dev, nlist_dev);

				// thermostat.apply_one(atom);
				// if (atom.bc) boundary.force_one(atom, neibs, nlist[i]);
			}
		);

		// convert back to Atoms structure
		thrust::copy(atoms_dev.begin(), atoms_dev.end(), atoms.begin());
	}

	void prepare(Atoms& atoms, Atoms& neibs, NeibsList& nlist) /*const*/
	{
		// convert to AtomsDev structure
		AtomsDev atoms_dev(atoms.begin(), atoms.end());

		DISPATCH_SOLVER_FOR_CUDA(solver_t, p_solver_dev);

		thrust::for_each(atoms_dev.begin(), atoms_dev.end(), [p_solver_dev] __device__ (auto&& atom)
		{
			atom.f  = Eigen::Matrix<double, 3, 1>::Zero();
			atom.V  = Eigen::Matrix<double, 3, 3>::Zero();
			atom.Ep = 0.0;

			// for EAM type
			atom.Em = 0.0;
			atom.nc = 0.0;

			// written assuming Single solver just for now
			// atom.rcut = potential()->get_rcut();
			atom.rcut = p_solver_dev->potential_impl()->get_rcut();

			atom.bc = 0x0;
			// boundary.prepare_one(atomit, neibs, nlist[i]);
		});

		// convert back to Atoms structure
		thrust::copy(atoms_dev.begin(), atoms_dev.end(), atoms.begin());
	}
#endif

	__host__ __device__ decltype(auto) potential() /*const*/
	{
		return static_cast<solver_t*>(this)->potential_impl();
	}
};

} // namespace mdcraft::solver::cuda
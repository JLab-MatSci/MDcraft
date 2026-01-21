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

	void prepare(Atoms& atoms, Atoms& neibs, NeibsList& nlist) /*const*/
	{

	}
};

template <typename S>
class TPairSolver : public PairSolver
{
public:
	TPairSolver() = default;

#ifdef mdcraft_ENABLE_CUDA_THRUST
	void forces(Atoms& atoms, Atoms& neibs, NeibsList& nlist) /*const*/
	{
		auto n_atoms = atoms.size();

		// convert to AtomsDev structure
		AtomsDev atoms_dev(atoms.begin(), atoms.end());

		auto first = thrust::make_zip_iterator(
			thrust::make_tuple(thrust::counting_iterator<std::size_t>(0), atoms_dev.begin())
		);
		auto last  = thrust::make_zip_iterator(
			thrust::make_tuple(thrust::counting_iterator<std::size_t>(n_atoms), atoms_dev.end())
		);

		// dispatch the correct solver version
		S* solver = static_cast<S*>(this);
		static_assert(std::is_trivially_copyable_v<std::remove_reference_t<S>>);
		static_assert(std::is_standard_layout_v<std::remove_reference_t<S>>);

		// convert into pointer accessible on device
		thrust::device_vector<S> v_solver_dev(1, *solver);
		S* p_solver_dev = thrust::raw_pointer_cast(v_solver_dev.data());

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
#endif

};

} // namespace mdcraft::solver::cuda
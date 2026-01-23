#pragma once

#include "base_potential.cuh"

namespace mdcraft { 
namespace solver { 
namespace potential {
namespace cuda {


class cuLJs : public TBasePotential<cuLJs> {
public:
#ifdef mdcraft_ENABLE_CUDA_THRUST
	__device__ void force_impl(
		data::Atom&             atom,
		data::Atom*             neibs,
		neibs::NeibsListOneDev& nlist
	) /*const*/
	{
		// std::vector<vector> rall(nlist.size());
		// std::vector<double> r2  (nlist.size());

		for (auto i = 0; i < nlist.size; ++i) 
		{
			const auto& j = nlist[i];
			auto& neib = neibs[i];

			vector const r_ji = neib.r - atom.r;
			double const r2ij = r_ji.squaredNorm();

			// r2[i]   = r2ij;
			// rall[i] = r_ji;
		}

		// std::size_t natoms = 0;
		// for (auto i = 0; i < nlist.size(); i++)
		// {
		// 	if (r2[i] > R2cut) continue;
		// 	if (r2[i] < 1e-15) continue;

		// 	r2[natoms] = r2[i];
		// 	rall[natoms] = rall[i];

		// 	++natoms;
		// }

		// atom.f += force(rall.data(), r2.data(), natoms);
	}
#endif

	constexpr Tag tag() const
	{
		return Tag::cuLJs;
	}

	double get_a() const
	{
		return a_;
	}

	void set_a(double a)
	{
		a_ = a;
		update_dev();
	}

	double get_sigma() const
	{
		return sigma_;
	}

	void set_sigma(double sigma)
	{
		sigma_ = sigma;
		update_dev();
	}

	double get_rcut() const
	{
		return rcut_;
	}

	void set_rcut(double rcut)
	{
		rcut_ = rcut;
		update_dev();
	}

protected:
	cuLJs(double a, double sigma, double rcut)
	  : a_(a)
	  , sigma_(sigma)
	  , rcut_(rcut)
	{}

	friend class TBasePotential<cuLJs>;

private:
	double a_, sigma_, rcut_;

	const double X2min_ = 1.2599210498948732;
};

} // cuda
} // potential
} // solver
} // mdcraft
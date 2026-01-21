#pragma once

#include "base_potential.cuh"

namespace mdcraft { 
namespace solver { 
namespace potential {
namespace cuda {

// pair Lennard-Jones potential V(r)[kJ/mol]
// X =  (r/rVr)^2 - X2min
// V(r) = 4* aVr * ( [rVr/r]^12 - [rVr/r]^6 + [aLJ3 + bLJ2*X]* X*X )

class cuLJs : public TBasePotential<cuLJs> {
public:
	cuLJs() {}

	// cuLJs(
	// 	double aVr,
	// 	double rVr,
	// 	double Rcutoff
	// );

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

private:
	
	double aVr, aVr2, rVr, rVr1, rVr2;
	double aLJ, aLJ3, bLJ, bLJ2;

	const double X2min = 1.2599210498948732;

};

} // cuda
} // potential
} // solver
} // mdcraft
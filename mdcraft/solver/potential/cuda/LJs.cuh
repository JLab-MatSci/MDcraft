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
		for (auto i = 0; i < nlist.size; ++i)
		{
			const auto& j = nlist[i];
			auto& neib = neibs[j];

			auto r_ij = dr::from_to(atom, neib);
			double const r2_ij = r_ij.squaredNorm();

			if (r2_ij < 1e-15) continue;

			if (r2_ij < rcut2_)
			{
				// V'(r) * (r_j - r_i) / r_ij: r_ij is cancelled in denominator here
				// and in numerator of V'(r^2)
				atom.f += dV_dr(r2_ij) * r_ij;
			}
		}
	}
#endif

	constexpr Tag tag() const
	{
		return Tag::cuLJs;
	}

	__host__ __device__ double get_a() const
	{
		return a_;
	}

	void set_a(double a)
	{
		a_ = a;
		update_dev();
	}

	__host__ __device__ double get_sigma() const
	{
		return sigma_;
	}

	void set_sigma(double sigma)
	{
		sigma_ = sigma;
		s2_    = sigma_*sigma_;
		s6_    = s2_*s2_*s2_;
		s12_   = s6_*s6_;
		update_coeffs_V();
		update_dev();
	}

	__host__ __device__ double get_rcut() const
	{
		return rcut_;
	}

	void set_rcut(double rcut)
	{
		rcut_  = rcut;
		rcut2_ = rcut_*rcut_;
		update_coeffs_V();
		update_dev();
	}

private:
	friend class TBasePotential<cuLJs>;

	cuLJs(double a, double sigma, double rcut)
	  : a_(a)
	  , sigma_(sigma)
	  , rcut_(rcut)
	  , rcut2_(rcut_*rcut_)
	{
		update_coeffs_V();
	}

	__host__ __device__ constexpr double X_V(double r2)
	{
		return r2 / s2_ - X2min_;
	}

	__host__ __device__ constexpr double sr2_V(double r)
	{
		auto src = sigma_ / rcut_;
		return src*src;
	}

	__host__ __device__ constexpr double sr6_V(double r)
	{
		return sr2_V(r)*sr2_V(r)*sr2_V(r);
	}

	__host__ __device__ constexpr double sr12_V(double r)
	{
		return sr6_V(r)*sr6_V(r);
	}

	/// Pair Lennard-Jones potential with corrections [kJ/mol]
	/// X = (r/sigma)^2 - X2min
	/// V(r) = 4 * a * ((sigma/r)^12 - (sigma/r)^6 + [b + c*X]*X*X)
	__host__ __device__ constexpr double V(double r)
	{
		auto sr = sigma_ / r;
		auto sr2 = sr*sr;
		auto sr6 = sr2*sr2*sr2;
		auto sr12 = sr6*sr6;
		auto X  = X_V(r*r);
		return 4.*a_ * (sr12 - sr6 + (b_ + c_*X)*X*X);  // ??? divide by two
	}

	__host__ __device__ constexpr double dV_dr(double r2)
	{
		auto r4  = r2*r2;
		auto r6  = r2*r4;
		auto r8  = r2*r6;
		auto r14 = r6*r8;
		auto X   = X_V(r2);
		return 4.*a_ * ( -12.*s12_/r14 + 6.*s6_/r8
					     + 2. * (2.*b_ + 3.*c_*X) * X / s2_);
	}

	void update_coeffs_V()
	{
		s2_    = sigma_*sigma_;
		s6_    = s2_*s2_*s2_;
		s12_   = s6_*s6_;

		Xc_ = X_V(rcut2_);
		auto sr6c = sr6_V(rcut_);
		auto D = (sr6c - 0.5) * Xc_ * s2_ / rcut2_;

		b_ = -3*(sr6c - 1. + 2.*D) * sr6c / (Xc_*Xc_);
		c_ =  2*(sr6c - 1. + 3.*D) * sr6c / (Xc_*Xc_*Xc_);
	}

private:
	double a_, sigma_, rcut_;

	double rcut2_, Xc_, b_, c_;

	double s2_, s6_, s12_;

	const double X2min_ = 1.2599210498948732;
};

} // cuda
} // potential
} // solver
} // mdcraft
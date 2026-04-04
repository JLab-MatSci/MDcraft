#include <mdcraft/solver/potential/LJs-cut.h>

namespace mdcraft {
namespace solver {
namespace potential {


LJsCut::LJsCut(
	double a,
	double s,
	double rcut
) : Base(rcut, ::mdcraft::solver::potential::type::pair)
  , a_(a)
  , s_(s)
  , s6_(std::pow(s, 6))
  , lj1_(2 * s6_)
  , lj2_(1.)
  , factor_lj_(48. * a_ * s6_)
{
}

void LJsCut::force(Atoms::iterator atomit, Atoms& neibs, NeibsListOne& nlist) {
	if (nlist.empty()) return;

	auto& atom = *atomit;

	auto& r_i = atom.r;
	double x{r_i(0)}, y{r_i(1)}, z{r_i(2)};

	auto sz = nlist.size();

	double fx{0.0}, fy{0.0}, fz{0.0};

	for (auto i = 0; i < sz; i++) {
		const auto& j = nlist[i];
		auto& neib = neibs[j];

		auto& r_j = neib.r;
		double xj{r_j(0)}, yj{r_j(1)}, zj{r_j(2)};

		auto drx = xj - x;
		auto dry = yj - y;
		auto drz = zj - z;
		auto r2ij = drx*drx + dry*dry + drz*drz;

		if ( (r2ij > R2cut) || (r2ij < 1e-15) ) continue;

		auto r2inv = 1.0 / r2ij;
		auto r6inv = r2inv * r2inv * r2inv;

		auto prefact = -0.5 * factor_lj_ * (lj1_ * r6inv - lj2_) * r6inv * r2inv /** 0.5*/;

		fx += drx * prefact;
		fy += dry * prefact;
		fz += drz * prefact;

		if (nlist.half())
		{
			neib.f(0) -= drx * prefact;
			neib.f(1) -= dry * prefact;
			neib.f(2) -= drz * prefact;
		}
	}

	atom.f(0) += fx;
	atom.f(1) += fy;
	atom.f(2) += fz;
}

void LJsCut::virial(Atoms::iterator atomit, Atoms& neibs, NeibsListOne& nlist) {
    if (nlist.empty()) return;
	auto& atom = *atomit;

	auto& r_i = atom.r;
	double x{r_i(0)}, y{r_i(1)}, z{r_i(2)};

	auto sz = nlist.size();

    double Vxx{0.0}, Vyy{0.0}, Vzz{0.0};
    double Vxy{0.0}, Vxz{0.0}, Vyz{0.0};

	for (auto i = 0; i < sz; i++) {
		const auto& j = nlist[i];
		auto& neib = neibs[j];

		auto& r_j = neib.r;
		double xj{r_j(0)}, yj{r_j(1)}, zj{r_j(2)};

		auto drx = xj - x;
		auto dry = yj - y;
		auto drz = zj - z;
		auto r2ij = drx*drx + dry*dry + drz*drz;

		if ( (r2ij > R2cut) || (r2ij < 1e-15) ) continue;

		auto r2inv = 1.0 / r2ij;
		auto r6inv = r2inv * r2inv * r2inv;

		auto prefact = -0.5 * factor_lj_ * (lj1_ * r6inv - lj2_) * r6inv * r2inv /** 0.5*/;

        auto fx = prefact * drx, fy = prefact * dry , fz = prefact * drz;

        Vxx += fx*drx;
        Vyy += fy*dry;
        Vzz += fz*drz;
        Vxy += fx*dry; // ???
        Vxz += fx*drz; // ???
        Vyz += fy*drz; // ???

		auto tmp = s6_ * r6inv;
		auto add_Ep = 0.5 * 4 * a_ * tmp * (tmp - 1);
		atom.Ep += add_Ep;

		if (nlist.half())
		{
			neib.V += 0.5 * matrix{
				{fx*drx, fx*dry, fx*drz},
				{fy*drx, fy*dry, fy*drz},
				{fz*drx, fz*dry, fz*drz}
			};

			neib.Ep += add_Ep;
		}
	}

    atom.V += 0.5 * matrix{
        {Vxx, Vxy, Vxz},
        {Vxy, Vyy, Vyz},
        {Vxz, Vyz, Vzz}
    };
}

} // namespace potential
} // namespace solver
} // namespace mdcraft

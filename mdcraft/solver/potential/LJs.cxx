#include <mdcraft/solver/potential/LJs.h>

namespace mdcraft { 
namespace solver { 
namespace potential {

LJs::LJs(
	double _aVr,
	double _rVr,
	double Rcutoff
) : Base(Rcutoff, ::mdcraft::solver::potential::type::pair)
{
	aVr  = _aVr;
	rVr  = _rVr;
	rVr2 =  rVr*rVr;
	rVr1 =  1.0/rVr2;
	aVr2 =  R2cut/rVr2 - X2min;
	aLJ3 =  (rVr2/R2cut)*(rVr2/R2cut)*(rVr2/R2cut);
	bLJ2 =  (aLJ3 - 0.5)*aVr2*rVr2/R2cut;
	aLJ  =  (aLJ3 - 1.0 + 2.0*bLJ2)*aLJ3/(aVr2*aVr2);
	bLJ  = -(aLJ3 - 1.0 + 3.0*bLJ2)*aLJ3/(aVr2*aVr2*aVr2);
	aLJ3 =  3.0*aLJ;
	bLJ2 =  2.0*bLJ;
	aVr2 =  -48.0*aVr/rVr2;
	aVr  =    2.0*aVr; // to get Epot/2
}
// potential U(r)
double LJs::value(const vector r) {
	double u, r2;
	r2 = r.squaredNorm();
	value(&r2, &u, 1ul);
	return u;
}
// virial U'(r)/r
matrix LJs::virial(const vector r) {
	double r2;
	matrix v;
	r2 = r.squaredNorm();
	virial(&r, &r2, &v, 1ul);
	return v;
}
// force
vector LJs::force(const vector r) {
	double r2 = r.squaredNorm();
	return force(&r, &r2, 1ul);
}

double LJs::value(const double* r2, std::size_t size) {
	double XX, FV, Ur, Uout = 0.0;

	for (std::size_t k = 0; k < size; ++k) {
		XX = rVr1 * r2[k] - X2min;
		FV = rVr2 / r2[k];
		Ur = FV * FV * FV;
		Uout += aVr  * ( Ur * (Ur - 1.0) - (aLJ3 + bLJ2*XX) * XX * XX); // V(r)/2
	}

	return Uout;
}

void LJs::value(const double* r2, double* Uout, std::size_t size) {
	double XX, FV, Ur;

	for (std::size_t k = 0; k < size; ++k) {
		XX = rVr1 * r2[k] - X2min;
		FV = rVr2 / r2[k];
		Ur = FV * FV * FV;
		Uout[k] = aVr  * ( Ur * (Ur - 1.0) - (aLJ3 + bLJ2*XX) * XX * XX); // V(r)/2
	}
}

void LJs::virial(const vector* r, const double* r2, matrix* Vout, std::size_t size) {
	double XX, FV, Ur;
	vector F;

	#pragma ivdep
	#pragma vector aligned
	for (std::size_t k = 0; k < size; ++k) {
        XX = rVr1 * r2[k] - X2min;
        FV = rVr2 / r2[k];
        Ur = FV * FV * FV;
        FV = aVr2 * (Ur * (Ur - 0.5)*FV + (aLJ + bLJ*XX)*XX); // V'(r)/r
        F  = FV * r[k];
        Vout[k] = 0.5 * F * r[k].transpose();
    }
    return;
}

void LJs::virial(Atoms::iterator atomit, Atoms& neibs, NeibsListOne& nlist) {
	if (nlist.half()) {
		auto& atom = *atomit;

		for (std::size_t i = 0; i < nlist.size(); ++i) {
			const auto j = nlist[i];
			auto& neib = neibs[j];

			const vector r_ji = neib.r - atom.r;
			const double r2ij = r_ji.squaredNorm();

			if (r2ij > R2cut || r2ij < 1e-15) continue;

			double XX = rVr1 * r2ij - X2min;
			double FV = rVr2 / r2ij;
			double Ur = FV * FV * FV;
			double FF = aVr2 * (Ur * (Ur - 0.5) * FV + (aLJ + bLJ * XX) * XX);
			const vector F = FF * r_ji;
			const matrix V = 0.5 * F * r_ji.transpose();
			const double E = aVr * (Ur * (Ur - 1.0) - (aLJ3 + bLJ2 * XX) * XX * XX);

			atom.V += V;
			atom.Ep += E;
			neib.V += V;
			neib.Ep += E;
		}

		return;
	}

	std::vector<vector> rall(nlist.size());
	std::vector<double> r2  (nlist.size());

	auto& atom = *atomit;

	matrix dummy_matrix;
	vector force_value;

	#pragma ivdep
	#pragma vector aligned
	for (auto i = 0; i < nlist.size(); i++) {
		auto& j = nlist[i];
		auto& neib = neibs[j];

		vector const r_ji = neib.r - atom.r;
		double const r2ij = r_ji.squaredNorm();

		r2[i]   = r2ij;
		rall[i] = r_ji;
	}

	std::size_t natoms = 0;
	for (auto i = 0; i < nlist.size(); i++){
		if (r2[i] > R2cut) continue;
		if (r2[i] < 1e-15) continue;

		r2[natoms] = r2[i];
		rall[natoms] = rall[i];
		++natoms;
	}

	atom.V  += virial(rall.data(), r2.data(), natoms);
	atom.Ep += value(r2.data(), natoms);
}

vector LJs::force(const vector* r, const double* r2, std::size_t size) {
	double XX, FV, Ur;
	vector Ftot = vector::Zero();

	#pragma ivdep
	#pragma vector aligned
	for (std::size_t k = 0; k < size; ++k) {
        XX = rVr1 * r2[k] - X2min;
        FV = rVr2 / r2[k];
        Ur = FV * FV * FV;
        FV = aVr2 * (Ur * (Ur - 0.5)*FV + (aLJ + bLJ*XX)*XX); // V'(r)/r
        Ftot += FV * r[k];
    }
    return Ftot;
}

matrix LJs::virial(const vector* r, const double* r2, std::size_t size) {
	double XX, FV, Ur;
	vector F = vector::Zero();
	matrix Vtot = matrix::Zero();

	#pragma ivdep
	#pragma vector aligned
	for (std::size_t k = 0; k < size; ++k) {
        XX = rVr1 * r2[k] - X2min;
        FV = rVr2 / r2[k];
        Ur = FV * FV * FV;
        FV = aVr2 * (Ur * (Ur - 0.5)*FV + (aLJ + bLJ*XX)*XX); // V'(r)/r
        F  = FV * r[k];
        Vtot += 0.5 * F * r[k].transpose();
    }
    return Vtot;
}

void LJs::force(const vector* r, vector* Fout, std::size_t size) {
	double XX, FV, Ur;
	for (std::size_t k = 0; k < size; ++k) {
		double r2 = r[k].squaredNorm();
        XX = rVr1 * r2 - X2min;
        FV = rVr2 / r2;
        Ur = FV * FV * FV;
        FV = aVr2 * (Ur * (Ur - 0.5)*FV + (aLJ + bLJ*XX)*XX); // V'(r)/r
        Fout[k] += FV * r[k];
    }
}

void LJs::force(Atoms::iterator atomit, Atoms& neibs, NeibsListOne& nlist) {
	if (nlist.half()) {
		auto& atom = *atomit;

		for (std::size_t i = 0; i < nlist.size(); ++i) {
			const auto j = nlist[i];
			auto& neib = neibs[j];

			const vector r_ji = neib.r - atom.r;
			const double r2ij = r_ji.squaredNorm();

			if (r2ij > R2cut || r2ij < 1e-15) continue;

			double XX = rVr1 * r2ij - X2min;
			double FV = rVr2 / r2ij;
			double Ur = FV * FV * FV;
			double FF = aVr2 * (Ur * (Ur - 0.5) * FV + (aLJ + bLJ * XX) * XX);
			const vector F = FF * r_ji;

			atom.f += F;
			neib.f -= F;
		}

		return;
	}

	std::vector<vector> rall(nlist.size());
	std::vector<double> r2  (nlist.size());

	auto& atom = *atomit;

	#pragma ivdep
	#pragma vector aligned
	for (auto i = 0; i < nlist.size(); i++) {
		const auto& j = nlist[i];
		auto& neib = neibs[j];

		vector const r_ji = neib.r - atom.r;
		double const r2ij = r_ji.squaredNorm();

		r2[i]   = r2ij;
		rall[i] = r_ji;
	}

	std::size_t natoms = 0;
	for (auto i = 0; i < nlist.size(); i++){
		if (r2[i] > R2cut) continue;
		if (r2[i] < 1e-15) continue;

		r2[natoms] = r2[i];
		rall[natoms] = rall[i];

		++natoms;
	}

	atom.f += force(rall.data(), r2.data(), natoms);
}

} // potential
} // solver
} // mdcraft

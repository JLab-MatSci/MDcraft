#pragma once
#include <vector>
#include <mutex>
#include <cstddef>

#include <mdcraft/solver/potential/base.h>

namespace mdcraft {
namespace solver {
namespace potential {

class EAMMOD1044 : public Base {
public:
    EAMMOD1044(
		const double       Rcutoff,
		const double       aVr,
		const double       bVr,
		const double       rVr,
		const std::size_t  arrsz,
		const double*      V,
		const double*      W,
		const double*      C,
		const double*      D,
		const double*      F,
		const double*      G
	);

	// potential U(r)
	double value(const double r);
	double value(const vector r);
	void   value(const double r2, double& Uout);

	// virial U'(r)/r
	matrix virial(const vector r, const double dF);
	void   virial(Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist);
	void   virial(const vector  r, const double  r2, const double  dF, matrix& Vout);

	// force
	vector force(const vector r, const double dF);
	void   force(Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist);
	void   force(const vector r, const double r2, const double dF, vector& Fout);

	// charge density
	void   density(Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist);
	void   density(const double r2, double& n);

	// embedding energy
	void embeddingEnergy(const double nr, double& EF);
	void embeddingEnergyDeriv(const double nr, double& dF);

	//temp
	void dRho(const double r2, double& dn);

private:
	double            R2cut;
	double            aVr, bVr, rVr;
	static const std::size_t arrsz = 32;
	double            V[arrsz],
	                  W[arrsz],
	                  C[arrsz],
	                  D[arrsz],
	                  F[arrsz],
	                  G[arrsz];

	mutable std::mutex m_;

	void   eFG(const double r2, double& F0, double& F1, double& G0, double& G1);
	void embGF(const double nr, double& F0, double& F1, double& G0, double& G1);

	void fill_relpos(Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist,
                     std::vector<vector>& rall, std::vector<double>& r2);
	void HydrogenCase(Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist,
                      const std::vector<vector>& rall, const std::vector<double>& r2);
	void HeavyAtomCase(Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist,
                      const std::vector<vector>& rall, const std::vector<double>& r2);

	//std::vector<double> Femb, Rho, Vpar, dFemb, dRho, dVpar;
	//int Nr1, Nrh1, N;
	// double dr, drh, dr1, drh1;
};

} // potential
} // solver
} // mdcraft
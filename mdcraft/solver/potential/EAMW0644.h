#pragma once
#include <vector>
#include <mutex>
#include <cstddef>
#include <vector>

#include <mdcraft/data/vector.h>
#include <mdcraft/solver/potential/base.h>

namespace mdcraft {
namespace solver {
namespace potential {

class EAMW0644 : public Base {
public:
	 EAMW0644(
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
	void   value(const double* r2, double* Uout, std::size_t size);
	double value(const double* r2, std::size_t size);

	// virial U'(r)/r
	matrix virial(const vector r, const double dF);
	matrix virial(const vector*  r, const double*  r2, const double*  dF, std::size_t    size);
	void   virial(const vector* r, const double* r2, const double* dF, matrix* Vout, std::size_t size);
	void   virial(Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist);

	// force
	vector force(const vector r, const double dF);
	vector force(const vector* r, const double* r2, const double* dF, std::size_t   size) ;
	void   force(const vector* r, const double* r2, const double* dF, vector* Fout, std::size_t size);
	void   force(Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist);

	// charge density
	double density(const double* r2, std::size_t size);
	void   density(Atoms::iterator self, Atoms& neibs, NeibsListOne& nlist);

	// embedding energy
	void   embeddingEnergy(const double*  nr, double* EF, double* dF,std::size_t size);

	//temp
	void dRho(const double* r2, double* dn, std::size_t size);

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

	//std::vector<double> Femb, Rho, Vpar, dFemb, dRho, dVpar;
	//int Nr1, Nrh1, N;
	// double dr, drh, dr1, drh1;
};

} // potential
} // solver
} // mdcraft
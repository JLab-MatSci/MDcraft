#pragma once
#include <mdcraft/solver/potential/base.h>

#include <chrono>

namespace mdcraft { 
namespace solver { 
namespace potential {

// pair Lennard-Jones potential V(r)[kJ/mol]
// V(r) = 4 * a * ( [s/r]^12 - [s/r]^6 )

class LJsCut : public Base {
public:
	LJsCut(
		double a,
		double s,
		double rcut
	);
	
	void   virial(Atoms::iterator atom, Atoms& neibs, NeibsListOne& nlist) override;
	void   force(Atoms::iterator atom, Atoms& neibs, NeibsListOne& nlist) override;

private:
	double a_, s_, s6_, lj1_, lj2_, factor_lj_;
};

}
}
}

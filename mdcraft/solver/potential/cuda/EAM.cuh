#pragma once

#include "base_potential.cuh"

namespace mdcraft { 
namespace solver { 
namespace potential {
namespace cuda {

// pair Lennard-Jones potential V(r)[kJ/mol]
// X =  (r/rVr)^2 - X2min
// V(r) = 4* aVr * ( [rVr/r]^12 - [rVr/r]^6 + [aLJ3 + bLJ2*X]* X*X )

class cuEAM : public TBasePotential<cuEAM> {
public:
	cuEAM() {}

	__device__ void forces(
	/* take info on all the atoms:
		total number,
		all positions,
		all types,
		neighbors' info   !	
	 */
	) /*const*/
	{
		printf("__device__ cuEAM::forces\n");
		// test use of kernel:
		// kernel_();

		/* 1. evaluate tid from CUDA */

		/* 2. no loop, but if (ii < n_total) */

		/* 3. evaluate energy, force, virial and save into some form on GPU */

		/* 4. move to CPU */
	}

	constexpr Tag tag() /*const*/
	{
		return Tag::cuLJs;
	}

private:

};


} // cuda
} // potential
} // solver
} // mdcraft
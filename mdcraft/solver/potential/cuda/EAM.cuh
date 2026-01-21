#pragma once

#include "base_potential.cuh"

namespace mdcraft { 
namespace solver { 
namespace potential {
namespace cuda {

class cuEAM : public TBasePotential<cuEAM> {
public:
	cuEAM() {}

#ifdef mdcraft_ENABLE_CUDA_THRUST
	__device__ void force_impl(
		data::Atom&             atom,
		data::Atom*             neibs, 
		neibs::NeibsListOneDev& nlist
	) /*const*/
	{
		
	}
#endif

	constexpr Tag tag() /*const*/
	{
		return Tag::cuEAM;
	}

private:

};


} // cuda
} // potential
} // solver
} // mdcraft
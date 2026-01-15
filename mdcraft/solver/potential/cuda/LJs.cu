#include "LJs.cuh"

namespace mdcraft { 
namespace solver { 
namespace potential {
namespace cuda {

// __device__ void cuLJs::forces(
// 	/* take info on all the atoms:
// 		total number,
// 		all positions,
// 		all types,
// 		neighbors' info   !	
// 	 */
// ) /*const*/
// {
// 	printf("__device__ cuLJs::forces\n");

// 	// test use of kernel:
// 	// kernel_();

// 	/* 1. evaluate tid from CUDA */

// 	/* 2. no loop, but if (ii < n_total) */

// 	/* 3. evaluate energy, force, virial and save into some form on GPU */

// 	/* 4. move to CPU */
// }

} // cuda
} // potential
} // solver
} // mdcraft
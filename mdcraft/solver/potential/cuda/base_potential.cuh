#pragma once

#include <string>
#include <functional>

namespace mdcraft { 
namespace solver { 
namespace potential {
namespace cuda {

// using kernel_type = __device__ void(*)();

// __device__ void test_kernel(char* a, std::string s)
// {

// }

// void(*pp)(char * a, std::string s) = test_kernel;

enum class Tag
{
	NonImpl,
	cuLJs
};

class BasePotential
{
public:
	constexpr Tag tag() const
	{
		return Tag::NonImpl;
	}
};

template <typename P>
class TBasePotential : public BasePotential
{
public:
	// BasePotential(kernel_type k)
	// 	: kernel_(k) {}

	__device__ void forces() /*const*/
	{
		static_cast<P*>(this)->forces();
	}

	constexpr Tag tag() /*const*/
	{
		return static_cast<P*>(this)->tag();
	}
protected:
	//Kernel kernel_;
};

template <typename P>
__global__ void k_forces(TBasePotential<P>* pot)
{
	pot->forces();
}

// class BasePotential : public TBasePotential<BasePotential>
// {

// };

// some tests...
/*
class DynPolyBase
{
public:
	virtual __device__ void forces() = 0;
};
*/

} // cuda
} // potential
} // solver
} // mdcraft
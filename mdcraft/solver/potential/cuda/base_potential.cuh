#pragma once

#include <mdcraft/configuration.h>

#include <mdcraft/data/atom.h>
#include <mdcraft/solver/stepper/base.h>

namespace mdcraft { 
namespace solver { 
namespace potential {
namespace cuda {

using NeibsList = neibs::List;
using Atoms     = mdcraft::data::Atoms;
using AtomsDev  = mdcraft::data::AtomsDev;
using NeibsListOne = ::mdcraft::neibs::ListOne;

enum class Tag
{
	NonImpl,
	cuLJs,
	cuEAM
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
#ifdef mdcraft_ENABLE_CUDA_THRUST
	__device__ void force(
		data::Atom&             atom,
		data::Atom*             neibs, 
		neibs::NeibsListOneDev& nlist
	) /*const*/
	{
		static_cast<P*>(this)->force_impl(atom, neibs, nlist);
	}
#endif

	constexpr Tag tag() /*const*/
	{
		return static_cast<P*>(this)->tag();
	}
};

} // cuda
} // potential
} // solver
} // mdcraft
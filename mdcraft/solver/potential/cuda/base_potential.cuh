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
	/// Manually free memory for both CPU and GPU pointers to P object.
	struct Deleter
	{
		void operator()(P* p) const
		{
			cudaFree(p->ptr_dev());
			delete p;
		}
	};

public:
	/// Custom wrapper for RAII support
	using PtrHolder = std::unique_ptr<P, Deleter>;

	/// Factory with RAII support both on CPU and GPU via custom PtrHolder.
	/// Each derived class of type P hides its ctor and friends with this type.
	template <typename... Args>
	static auto create(Args&&... args)
	{
		static_assert(std::is_trivially_copyable_v<std::remove_reference_t<P>>);

		// call ctor of the specific potential
		PtrHolder pot(new P{std::forward<Args>(args)...}, Deleter{});

		// allocate pointer-to-device-memory for created pot
		cudaMalloc(&pot->potential_dev_, sizeof(P));

		// copy to device
		cudaMemcpy(pot->potential_dev_, pot.get(), sizeof(P), cudaMemcpyHostToDevice);

		assert(pot->ptr_dev() != nullptr);

		return pot;
	}

	/// NO DTOR TO FREE MEMORY. ELSE WE LOOSE TRIVIALLY COPYABLE CLASS
	/// See Deleter struct below
	/*~TSingle()
	{
		cudaFree(potential_dev_);
	}*/

	P* ptr_dev() const
	{
		return potential_dev_;
	}

	/// Used in all setter methods of derived classes.
	void update_dev()
	{
		cudaMemcpy(potential_dev_, static_cast<P*>(this), sizeof(P), cudaMemcpyHostToDevice);
	}

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

protected:
	/// Simple way to additionally save device representation of P object.
	P* potential_dev_ = nullptr;
};

} // cuda
} // potential
} // solver
} // mdcraft
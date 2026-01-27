// #pragma once

// #include <mdcraft/external/variant/variant.h>

// #include "pair.h"
// #include <mdcraft/solver/potential/cuda/LJs.cuh>

// namespace mdcraft::solver::cuda {

// using Boundary  = boundary::Base;
// using cuLJs     = mdcraft::solver::potential::cuda::cuLJs;

// class Single : public TPairSolver<Single>
// {
// public:
// 	Single(
// 		Boundary&   boundary   = boundary::dummy_boundary,
// 		Thermostat& thermostat = thermostat::dummy_thermostat
// 	);

// 	void forces(Atoms& atoms, Atoms& neibs, NeibsList& nlist) /*const*/
// 	{
// 		pot_visitor visitor{atoms, neibs, nlist};
// 		variant::apply_visitor(visitor, potential_);
// 		//variant::apply_visitor([] (const auto& p) { p.forces(); }, potential);
// 	}

// private:
// 	variant::variant<cuLJs> potential_;
// 	/**<
// 		\~russian Межатомный потенциал.
// 		\~english Interatomic potential.
// 		\~
// 	*/

// 	struct pot_visitor
// 	{
// 		__device__ void operator()(const cuLJs& pot)
// 		{
// 			return pot.forces(/*atoms, neibs, nlist*/);
// 		}

// 		Atoms& atoms;
// 		Atoms& neibs;
// 		NeibsList& nlist;
// 	};
// };

// } // namespace mdcraft::solver::cuda

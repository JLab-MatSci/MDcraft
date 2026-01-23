#pragma once

#include <mdcraft/solver/cuda/isolver.h>

#ifdef mdcraft_ENABLE_CUDA_VARIANT
	#include <mdcraft/solver/cuda/single_var.h>

#elif defined(mdcraft_ENABLE_CUDA_CRTP)
	#include <mdcraft/solver/cuda/single_crtp.h>
	#include <mdcraft/solver/cuda/multi_crtp.h>

	#include <mdcraft/solver/potential/cuda/LJs.cuh>
	#include <mdcraft/solver/potential/cuda/EAM.cuh>

namespace {

// tag
	using Tag = mdcraft::solver::potential::cuda::Tag;

// solvers
	using BaseSolver = mdcraft::solver::cuda::BaseSolver;
	using PairSolverBase = mdcraft::solver::cuda::PairSolver;
	using MultiSolver  = mdcraft::solver::cuda::Multi;

	template <typename S>
	using PairSolver = mdcraft::solver::cuda::TPairSolver<S>;

	template <typename P>
	using SingleCUDA = mdcraft::solver::cuda::TSingle<P>;

	template <typename P>
	using MultiCUDA = mdcraft::solver::cuda::TMulti<P>;

// potentials
	using cuLJs = mdcraft::solver::potential::cuda::cuLJs;
	using cuEAM = mdcraft::solver::potential::cuda::cuEAM;

	using BasePotential = mdcraft::solver::potential::cuda::BasePotential;

	template <typename P>
	using TBasePotential = mdcraft::solver::potential::cuda::TBasePotential<P>;
#endif

}
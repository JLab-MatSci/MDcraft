#pragma once

#ifdef mdcraft_ENABLE_CUDA_VARIANT
	#include <mdcraft/solver/cuda/single_var.h>

#elif defined(mdcraft_ENABLE_CUDA_CRTP)
	#include <mdcraft/solver/cuda/single_crtp.h>
	#include <mdcraft/solver/cuda/multi_crtp.h>

	#include <mdcraft/solver/potential/cuda/LJs.cuh>
	#include <mdcraft/solver/potential/cuda/EAM.cuh>

// tag
	using Tag = mdcraft::solver::potential::cuda::Tag;

// solvers
	using BaseSolver = mdcraft::solver::cuda::BaseSolver;
	using SingleSolver = mdcraft::solver::cuda::Single;
	using MultiSolver  = mdcraft::solver::cuda::Multi;

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
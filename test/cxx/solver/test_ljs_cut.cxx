#include <gtest/gtest.h>

#include <cmath>

#include <mdcraft/solver/potential/LJs-cut.h>
#include <mdcraft/solver/single.h>
#include <mdcraft/solver/boundary/base.h>
#include <mdcraft/tools/threads.h>

namespace mdcraft::test {

namespace {

double lj_pair_energy_half(double epsilon, double sigma, double distance) {
  const double sr = sigma / distance;
  const double sr6 = std::pow(sr, 6);
  return 0.5 * 4.0 * epsilon * sr6 * (sr6 - 1.0);
}

double lj_force_prefactor(double epsilon, double sigma, double distance) {
  const double sigma6 = std::pow(sigma, 6);
  const double factor = 48.0 * epsilon * sigma6;
  const double r2 = distance * distance;
  const double r2inv = 1.0 / r2;
  const double r6inv = r2inv * r2inv * r2inv;
  return -0.5 * factor * ((2.0 * sigma6) * r6inv - 1.0) * r6inv * r2inv;
}

}  // namespace

TEST(LJsCutTest, SolverPrepareResetsStateAndSetsCutoff) {
  tools::Threads threads(1);
  solver::potential::LJsCut potential(1.03120074782442750, 0.33841043857528683, 0.8125);
  solver::Single solver(potential, solver::boundary::dummy_boundary, solver::thermostat::dummy_thermostat, threads);

  data::Atoms atoms(1);
  atoms[0].m = 40.0;
  atoms[0].f << 1.0, 2.0, 3.0;
  atoms[0].Ep = 5.0;
  atoms[0].Em = 7.0;
  atoms[0].nc = 9.0;
  atoms[0].V = data::matrix::Identity();

  neibs::List nlist(1);
  solver.prepare(atoms, atoms, nlist);

  EXPECT_EQ(solver.stages.size(), 1U);
  EXPECT_DOUBLE_EQ(atoms[0].rcut, potential.rcut());
  EXPECT_TRUE(atoms[0].f.isZero(0.0));
  EXPECT_TRUE(atoms[0].V.isZero(0.0));
  EXPECT_DOUBLE_EQ(atoms[0].Ep, 0.0);
}

TEST(LJsCutTest, ForceAndVirialMatchAnalyticalPairValues) {
  constexpr double epsilon = 1.03120074782442750;
  constexpr double sigma = 0.33841043857528683;
  constexpr double rcut = 0.8125;
  constexpr double distance = 0.40;

  solver::potential::LJsCut potential(epsilon, sigma, rcut);

  data::Atoms atoms(2);
  atoms[0].r << 0.0, 0.0, 0.0;
  atoms[1].r << distance, 0.0, 0.0;

  neibs::ListOne nlist0{1U};
  neibs::ListOne nlist1{0U};

  potential.force(atoms.begin(), atoms, nlist0);
  potential.force(atoms.begin() + 1, atoms, nlist1);
  potential.virial(atoms.begin(), atoms, nlist0);
  potential.virial(atoms.begin() + 1, atoms, nlist1);

  const double expected_prefactor = lj_force_prefactor(epsilon, sigma, distance);
  const double expected_fx = expected_prefactor * distance;
  const double expected_energy = 2.0 * lj_pair_energy_half(epsilon, sigma, distance);
  const double expected_virial = expected_fx * distance;

  EXPECT_NEAR(atoms[0].f.x(), expected_fx, 1e-12);
  EXPECT_NEAR(atoms[1].f.x(), -expected_fx, 1e-12);
  EXPECT_NEAR(atoms[0].f.y(), 0.0, 1e-12);
  EXPECT_NEAR(atoms[1].f.y(), 0.0, 1e-12);
  EXPECT_NEAR(atoms[0].Ep + atoms[1].Ep, expected_energy, 1e-12);
  EXPECT_NEAR(atoms[0].V.trace() + atoms[1].V.trace(), expected_virial, 1e-12);
}

TEST(LJsCutTest, HalfListUpdatesBothAtomsOnce) {
  constexpr double epsilon = 1.03120074782442750;
  constexpr double sigma = 0.33841043857528683;
  constexpr double rcut = 0.8125;
  constexpr double distance = 0.40;

  solver::potential::LJsCut potential(epsilon, sigma, rcut);

  data::Atoms atoms(2);
  atoms[0].r << 0.0, 0.0, 0.0;
  atoms[1].r << distance, 0.0, 0.0;

  neibs::ListOne nlist0(true);
  nlist0.push_back(0U);
  nlist0.push_back(1U);
  neibs::ListOne empty_half(true);

  potential.force(atoms.begin(), atoms, nlist0);
  potential.force(atoms.begin() + 1, atoms, empty_half);
  potential.virial(atoms.begin(), atoms, nlist0);
  potential.virial(atoms.begin() + 1, atoms, empty_half);

  const double expected_prefactor = lj_force_prefactor(epsilon, sigma, distance);
  const double expected_fx = expected_prefactor * distance;
  const double expected_energy = 2.0 * lj_pair_energy_half(epsilon, sigma, distance);
  const double expected_virial = expected_fx * distance;

  EXPECT_NEAR(atoms[0].f.x(), expected_fx, 1e-12);
  EXPECT_NEAR(atoms[1].f.x(), -expected_fx, 1e-12);
  EXPECT_NEAR(atoms[0].Ep + atoms[1].Ep, expected_energy, 1e-12);
  EXPECT_NEAR(atoms[0].V.trace() + atoms[1].V.trace(), expected_virial, 1e-12);
}

}  // namespace mdcraft::test

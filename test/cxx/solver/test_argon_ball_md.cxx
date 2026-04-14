#include <gtest/gtest.h>

#include <memory>

#include <mdcraft/solver/boundary/base.h>
#include <mdcraft/solver/potential/LJs-cut.h>
#include <mdcraft/solver/single.h>
#include <mdcraft/solver/stepper/verlet.h>
#include <mdcraft/tools/threads.h>

#include <test/cxx/support/md_test_utils.h>

namespace mdcraft::test {

class ArgonBallMdTest : public ::testing::Test {
protected:
  static constexpr double kArgonMass = 40.0;
  static constexpr double kTimeStep = 4e-3;
  static constexpr double kNeighborSkinFactor = 10.125 / 8.125;

  tools::Threads threads{1};
  solver::potential::LJsCut potential{1.03120074782442750, 0.33841043857528683, 0.8125};
  solver::Single solver{
    potential,
    solver::boundary::dummy_boundary,
    solver::thermostat::dummy_thermostat,
    threads
  };
  solver::stepper::Verlet stepper{solver, threads};
  Atoms atoms;
  Domain domain;
  std::unique_ptr<VerletList> nlist;

  void SetUp() override {
    auto data = load_lammps_atomic_data(argon_ball_path());
    atoms = std::move(data.atoms);
    domain = data.domain;
    set_common_atom_properties(atoms, kArgonMass, potential.rcut(), kNeighborSkinFactor);
    zero_velocities(atoms);
    nlist = std::make_unique<VerletList>(atoms, atoms, domain, threads);
  }
};

TEST_F(ArgonBallMdTest, LoadsExpectedAtomicSystemFromLammpsData) {
  ASSERT_EQ(atoms.size(), 130124U);
  EXPECT_NEAR(domain.xsize(), 31.5, 1e-12);
  EXPECT_NEAR(domain.ysize(), 31.5, 1e-12);
  EXPECT_NEAR(domain.zsize(), 31.5, 1e-12);
  EXPECT_NEAR(atoms.front().m, 40.0, 1e-12);
  EXPECT_EQ(atoms.front().uid, 1U);
  EXPECT_EQ(atoms.back().uid, 130124U);
}

TEST_F(ArgonBallMdTest, ShortMdRunMatchesReferencePotentialAndVirialOnStepTwo) {
  const auto history = run_md_steps(solver, stepper, atoms, *nlist, kTimeStep, 2, 1);

  ASSERT_EQ(history.size(), 3U);
  EXPECT_NEAR(history[2].potential_energy, -1.0244383150033894e+06, 1.0e-3);
  EXPECT_NEAR(history[2].virial, -1.1164231473675424e+05, 1.0e-3);
  EXPECT_NEAR(history[2].total_energy, history[0].total_energy, 2.0);
}

TEST_F(ArgonBallMdTest, ShortMdRunProducesStableParticleLevelRegressionValues) {
  const auto history = run_md_steps(solver, stepper, atoms, *nlist, kTimeStep, 2, 1);
  (void)history;

  const std::size_t atom1 = find_atom_index_by_uid(atoms, 1U);
  const std::size_t atom2 = find_atom_index_by_uid(atoms, 2U);

  EXPECT_NEAR(atoms[atom1].r.x(), -3.9184367989026225, 1.0e-12);
  EXPECT_NEAR(atoms[atom1].r.y(), 2.6440651643797461, 1.0e-12);
  EXPECT_NEAR(atoms[atom1].r.z(), 5.2690712254285517, 1.0e-12);
  EXPECT_NEAR(atoms[atom1].v.z(), 5.79013209455442e-04, 1.0e-15);
  EXPECT_NEAR(atoms[atom1].f.z(), 2.8929546164515996, 1.0e-12);
  EXPECT_NEAR(atoms[atom1].Ep, -4.0629298955784474, 1.0e-12);
  EXPECT_NEAR(atoms[atom1].V.trace(), -3.3178872346261723e-01, 1.0e-12);

  EXPECT_NEAR(atoms[atom2].r.z(), 5.7940689066038047, 1.0e-12);
  EXPECT_NEAR(atoms[atom2].v.x(), -1.4268458571713898e-03, 1.0e-15);
  EXPECT_NEAR(atoms[atom2].f.x(), -7.1268568973069453, 1.0e-12);
  EXPECT_NEAR(atoms[atom2].Ep, -4.2002585997144832, 1.0e-12);
  EXPECT_NEAR(atoms[atom2].V.trace(), 4.2842826328841377e-01, 1.0e-12);
}

}  // namespace mdcraft::test

#include <gtest/gtest.h>

#include <memory>

#include <mdcraft/neibs/verlet-list.h>
#include <mdcraft/solver/boundary/periodic.h>
#include <mdcraft/solver/potential/EAM.h>
#include <mdcraft/solver/single.h>
#include <mdcraft/solver/stepper/verlet.h>
#include <mdcraft/tools/threads.h>

#include <test/cxx/support/md_test_utils.h>

namespace mdcraft::test {

class AluminumBulkEamTest : public ::testing::Test {
protected:
  static constexpr double kAlMass = 26.9815384;
  static constexpr double kLatticeConstant = 0.4032;
  static constexpr double kTimeStep = 4e-3;
  static constexpr double kNeighborSkinFactor = 1.25;

  tools::Threads threads{1};
  solver::potential::EAM potential{aluminum_eam_path().string()};
  Atoms atoms;
  Domain domain;
  std::unique_ptr<solver::boundary::Periodic> periodic;
  std::unique_ptr<solver::Single> solver_ptr;
  std::unique_ptr<solver::stepper::Verlet> stepper_ptr;
  std::unique_ptr<neibs::VerletList> nlist;

  void SetUp() override {
    auto [bulk_atoms, bulk_domain] = make_fcc_bulk(3, 3, 3, kLatticeConstant, kAlMass);
    atoms = std::move(bulk_atoms);
    domain = bulk_domain;
    set_common_atom_properties(atoms, kAlMass, potential.rcut(), kNeighborSkinFactor);
    assign_deterministic_velocities(atoms, 5.0e-3);

    periodic = std::make_unique<solver::boundary::Periodic>(domain);
    solver_ptr = std::make_unique<solver::Single>(
      potential,
      *periodic,
      solver::thermostat::dummy_thermostat,
      threads
    );
    stepper_ptr = std::make_unique<solver::stepper::Verlet>(*solver_ptr, threads);
    nlist = std::make_unique<neibs::VerletList>(atoms, atoms, domain, threads);
  }
};

TEST_F(AluminumBulkEamTest, EamSolverUsesTwoStagesAndPeriodicNeighbors) {
  nlist->update(-1.0);

  ASSERT_EQ(solver_ptr->stages.size(), 2U);
  EXPECT_EQ(atoms.size(), 108U);
  EXPECT_FALSE(nlist->get().empty());
  EXPECT_GT(nlist->operator[](0).size(), 12U);
}

TEST_F(AluminumBulkEamTest, ShortPeriodicBulkRunConservesEnergyWithoutThermostat) {
  const auto history = run_md_steps(*solver_ptr, *stepper_ptr, atoms, *nlist, kTimeStep, 5, 1);

  ASSERT_EQ(history.size(), 6U);
  EXPECT_NEAR(history.back().total_energy, history.front().total_energy, 40.0);
  EXPECT_NEAR(history.back().com_velocity.norm(), 0.0, 1.0e-12);
  EXPECT_TRUE(std::isfinite(history.back().potential_energy));
  EXPECT_TRUE(std::isfinite(history.back().embedding_energy));
  EXPECT_TRUE(std::isfinite(history.back().virial));
}

TEST_F(AluminumBulkEamTest, ShortPeriodicBulkRunKeepsRepresentativeAtomsFiniteAndPhysical) {
  const auto history = run_md_steps(*solver_ptr, *stepper_ptr, atoms, *nlist, kTimeStep, 2, 1);

  EXPECT_NEAR(history.back().total_energy, -3.4489006959225437e+04, 1.0e-6);
  EXPECT_NEAR(history.back().potential_energy, -6.7190981865001777e+03, 1.0e-6);
  EXPECT_NEAR(history.back().embedding_energy, -2.7770493824864527e+04, 1.0e-6);
  EXPECT_NEAR(history.back().kinetic_energy, 5.8505213926282473e-01, 1.0e-9);
  EXPECT_NEAR(history.back().virial, 7.4272691119381875e+01, 1.0e-9);

  EXPECT_NEAR(atoms[0].r.x(), 6.0469801558679681e-01, 1.0e-12);
  EXPECT_NEAR(atoms[0].r.y(), 6.0474258946107018e-01, 1.0e-12);
  EXPECT_NEAR(atoms[0].r.z(), 6.0477602485847293e-01, 1.0e-12);
  EXPECT_NEAR(atoms[0].v.x(), -1.3221213313405629e-02, 1.0e-15);
  EXPECT_NEAR(atoms[0].f.x(), -1.1358391912725722e+01, 1.0e-12);
  EXPECT_NEAR(atoms[0].Ep, -6.3230321656455359e+01, 1.0e-9);
  EXPECT_NEAR(atoms[0].Em, -2.5721082474679707e+02, 1.0e-9);
  EXPECT_NEAR(atoms[0].V.trace(), 1.3557390418206403e+01, 1.0e-9);
}

}  // namespace mdcraft::test

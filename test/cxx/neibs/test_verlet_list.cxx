#include <gtest/gtest.h>

#include <algorithm>
#include <stdexcept>

#include <mdcraft/neibs/grid.h>
#include <mdcraft/neibs/verlet-list.h>
#include <mdcraft/tools/threads.h>

namespace mdcraft::test {

namespace {

data::Atom make_atom(double x, double y, double z, double rcut, double rns) {
  data::Atom atom{};
  atom.r << x, y, z;
  atom.rcut = rcut;
  atom.rns = rns;
  atom.m = 1.0;
  return atom;
}

}  // namespace

TEST(VerletListTest, FindsNeighborsWithoutPeriodicity) {
  tools::Threads threads(1);

  data::Atoms atoms(3);
  atoms[0].r << 0.20, 0.0, 0.0;
  atoms[1].r << 0.85, 0.0, 0.0;
  atoms[2].r << 2.50, 0.0, 0.0;

  for (auto& atom : atoms) {
    atom.rns = 0.80;
    atom.rcut = 0.80;
    atom.m = 1.0;
  }

  lattice::Domain domain(0.0, 3.0, -1.0, 1.0, -1.0, 1.0);
  neibs::VerletList nlist(atoms, atoms, domain, threads);
  nlist.update(-1.0);

  auto first = nlist[0];
  std::sort(first.begin(), first.end());

  EXPECT_EQ(nlist[0].size(), 2U);
  EXPECT_EQ(nlist[1].size(), 2U);
  EXPECT_EQ(nlist[2].size(), 1U);
  EXPECT_EQ(first[0], 0U);
  EXPECT_EQ(first[1], 1U);
  EXPECT_EQ(nlist[2][0], 2U);
}

TEST(VerletListTest, FindsNeighborsThroughPeriodicBoundary) {
  tools::Threads threads(1);

  data::Atoms atoms(2);
  atoms[0].r << -0.48, 0.0, 0.0;
  atoms[1].r << 0.48, 0.0, 0.0;

  for (auto& atom : atoms) {
    atom.rns = 0.10;
    atom.rcut = 0.10;
    atom.m = 1.0;
  }

  lattice::Domain domain(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5);
  domain.set_periodic(0, true);

  neibs::VerletList nlist(atoms, atoms, domain, threads);
  nlist.update(-1.0);

  ASSERT_EQ(nlist[0].size(), 2U);
  ASSERT_EQ(nlist[1].size(), 2U);
  auto first = nlist[0];
  auto second = nlist[1];
  std::sort(first.begin(), first.end());
  std::sort(second.begin(), second.end());
  EXPECT_EQ(first[0], 0U);
  EXPECT_EQ(first[1], 1U);
  EXPECT_EQ(second[0], 0U);
  EXPECT_EQ(second[1], 1U);
}

TEST(VerletListTest, ListForReturnsNeighborsForQueryPoint) {
  tools::Threads threads(1);

  data::Atoms atoms;
  atoms.push_back(make_atom(0.15, 0.10, 0.10, 0.40, 0.40));
  atoms.push_back(make_atom(0.55, 0.10, 0.10, 0.40, 0.40));
  atoms.push_back(make_atom(1.60, 0.10, 0.10, 0.40, 0.40));

  lattice::Domain domain(0.0, 2.0, 0.0, 1.0, 0.0, 1.0);
  neibs::VerletList nlist(atoms, atoms, domain, threads);
  nlist.update(-1.0);

  data::vector query;
  query << 0.40, 0.10, 0.10;
  auto found = nlist.list_for(query);
  std::sort(found.begin(), found.end());

  ASSERT_EQ(found.size(), 2U);
  EXPECT_EQ(found[0], 0U);
  EXPECT_EQ(found[1], 1U);
}

TEST(VerletListTest, SortOrdersAtomsByOwningCell) {
  tools::Threads threads(1);

  data::Atoms atoms;
  atoms.push_back(make_atom(1.60, 0.10, 0.10, 0.40, 0.40));
  atoms.push_back(make_atom(0.15, 0.10, 0.10, 0.40, 0.40));
  atoms.push_back(make_atom(0.95, 0.10, 0.10, 0.40, 0.40));

  lattice::Domain domain(0.0, 2.0, 0.0, 1.0, 0.0, 1.0);
  neibs::VerletList nlist(atoms, atoms, domain, threads);

  auto sorted = nlist.sort(atoms);

  ASSERT_EQ(sorted.size(), 3U);
  EXPECT_LT(sorted[0].r.x(), sorted[1].r.x());
  EXPECT_LT(sorted[1].r.x(), sorted[2].r.x());
}

TEST(VerletListTest, UpdateWithBufferFactorUsesActualNeighborIndices) {
  tools::Threads threads(1);

  data::Atoms atoms;
  atoms.push_back(make_atom(0.10, 0.10, 0.10, 0.20, 0.35));
  atoms.push_back(make_atom(0.35, 0.10, 0.10, 0.60, 0.35));
  atoms.push_back(make_atom(1.80, 0.10, 0.10, 0.10, 0.35));

  lattice::Domain domain(0.0, 2.0, 0.0, 1.0, 0.0, 1.0);
  neibs::VerletList nlist(atoms, atoms, domain, threads);

  nlist.update(-1.0);
  nlist.update(2.0);

  EXPECT_NEAR(atoms[0].rns, 1.60, 1e-12);
  EXPECT_NEAR(atoms[1].rns, 2.40, 1e-12);
  EXPECT_NEAR(atoms[2].rns, 0.40, 1e-12);
}

TEST(VerletListTest, HalfListStoresEachLocalPairOnce) {
  tools::Threads threads(1);

  data::Atoms atoms(3);
  atoms[0] = make_atom(0.10, 0.10, 0.10, 0.35, 0.35);
  atoms[1] = make_atom(0.32, 0.10, 0.10, 0.35, 0.35);
  atoms[2] = make_atom(0.95, 0.10, 0.10, 0.35, 0.35);

  lattice::Domain domain(0.0, 2.0, 0.0, 1.0, 0.0, 1.0);
  neibs::VerletList nlist(atoms, atoms, domain, threads, true);
  nlist.update(-1.0);

  auto first = nlist[0];
  std::sort(first.begin(), first.end());

  ASSERT_EQ(nlist[0].size(), 2U);
  ASSERT_EQ(nlist[1].size(), 1U);
  ASSERT_EQ(nlist[2].size(), 1U);
  EXPECT_EQ(first[0], 0U);
  EXPECT_EQ(first[1], 1U);
  EXPECT_EQ(nlist[1][0], 1U);
  EXPECT_EQ(nlist[2][0], 2U);
}

TEST(VerletListTest, DecideUsesHalfSkinDisplacementCriterion) {
  tools::Threads threads(1);

  data::Atoms atoms(2);
  atoms[0] = make_atom(0.10, 0.10, 0.10, 0.60, 0.80);
  atoms[1] = make_atom(0.50, 0.10, 0.10, 0.60, 0.80);
  atoms[0].uid = 1;
  atoms[1].uid = 2;

  lattice::Domain domain(0.0, 2.0, 0.0, 1.0, 0.0, 1.0);
  neibs::VerletList nlist(atoms, atoms, domain, threads, true);
  nlist.update(-1.0);

  atoms[0].r.x() += 0.09;
  EXPECT_FALSE(nlist.decide());

  atoms[0].r.x() += 0.02;
  EXPECT_TRUE(nlist.decide());
}

TEST(VerletListTest, HalfListRejectsActiveThreadPool) {
  tools::Threads threads(2);
  data::Atoms atoms(1);
  atoms[0] = make_atom(0.10, 0.10, 0.10, 0.20, 0.30);

  lattice::Domain domain(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);

  EXPECT_THROW(neibs::VerletList(atoms, atoms, domain, threads, true), std::runtime_error);
}

}  // namespace mdcraft::test

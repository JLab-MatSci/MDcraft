#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <vector>

#include <mdcraft/neibs/grid.h>
#include <mdcraft/tools/threads.h>

namespace mdcraft::test {

namespace {

data::Atom make_atom(double x, double y, double z, double rns) {
  data::Atom atom{};
  atom.r << x, y, z;
  atom.rns = rns;
  atom.rcut = rns;
  atom.m = 1.0;
  return atom;
}

std::vector<neibs::cell_id_t> compact_cells(const std::array<neibs::cell_id_t, 27>& cells) {
  std::vector<neibs::cell_id_t> result;
  for (auto cell : cells) {
    if (cell == neibs::empty_cell) {
      break;
    }
    result.push_back(cell);
  }
  std::sort(result.begin(), result.end());
  result.erase(std::unique(result.begin(), result.end()), result.end());
  return result;
}

}  // namespace

TEST(GridTest, BuildCreatesExpectedLocalSubgridAndIndices) {
  tools::Threads threads(1);
  lattice::Domain domain(0.0, 4.0, 0.0, 4.0, 0.0, 4.0);

  data::Atoms atoms;
  atoms.push_back(make_atom(0.2, 0.2, 0.2, 1.0));
  atoms.push_back(make_atom(1.8, 1.1, 0.7, 1.0));
  atoms.push_back(make_atom(2.6, 1.9, 0.8, 1.0));

  neibs::Grid grid(domain, threads);
  grid.build(atoms);

  EXPECT_NEAR(grid.cell_size(), 1.0, 1e-12);
  EXPECT_EQ(grid.cells_number(), 6);

  const auto idx0 = grid.local_index3D(atoms[0].r);
  const auto idx1 = grid.local_index3D(atoms[1].r);
  const auto idx2 = grid.local_index3D(atoms[2].r);

  EXPECT_EQ(idx0[0], 0);
  EXPECT_EQ(idx0[1], 0);
  EXPECT_EQ(idx0[2], 0);
  EXPECT_EQ(idx1[0], 1);
  EXPECT_EQ(idx1[1], 1);
  EXPECT_EQ(idx1[2], 0);
  EXPECT_EQ(idx2[0], 2);
  EXPECT_EQ(idx2[1], 1);
  EXPECT_EQ(idx2[2], 0);

  EXPECT_EQ(grid.cell_index(idx0), 0);
  EXPECT_EQ(grid.cell_index(idx1), 4);
  EXPECT_EQ(grid.cell_index(idx2), 5);
}

TEST(GridTest, AdjacentCellsIncludeWrappedPeriodicNeighbors) {
  tools::Threads threads(1);
  lattice::Domain domain(-2.0, 2.0, -2.0, 2.0, -2.0, 2.0);
  domain.set_periodic(0, true);

  data::Atoms atoms;
  atoms.push_back(make_atom(-1.8, 0.0, 0.0, 1.0));
  atoms.push_back(make_atom(-0.2, 0.0, 0.0, 1.0));
  atoms.push_back(make_atom(0.2, 0.0, 0.0, 1.0));
  atoms.push_back(make_atom(1.8, 0.0, 0.0, 1.0));

  neibs::Grid grid(domain, threads);
  grid.build(atoms);

  const auto [cells, is_edge] = grid.adjacent_cells(atoms.front().r);
  const auto compact = compact_cells(cells);

  EXPECT_TRUE(is_edge[0]);
  EXPECT_EQ(compact.size(), 3U);
  EXPECT_EQ(compact[0], 0);
  EXPECT_EQ(compact[1], 1);
  EXPECT_EQ(compact[2], 3);
}

}  // namespace mdcraft::test

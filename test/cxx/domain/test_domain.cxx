#include <gtest/gtest.h>

#include <limits>

#include <mdcraft/lattice/domain.h>

namespace mdcraft::test {

TEST(DomainTest, BelongFitInPeriodAndShortestRespectPeriodicity) {
  lattice::Domain domain(-1.0, 1.0, -2.0, 2.0, -3.0, 3.0);
  domain.set_periodic(0, true);
  domain.set_periodic(1, true);

  data::vector point(1.25, -2.25, 0.5);
  domain.fit_in_period(point);

  EXPECT_NEAR(point.x(), -0.75, 1e-12);
  EXPECT_NEAR(point.y(), 1.75, 1e-12);
  EXPECT_NEAR(point.z(), 0.5, 1e-12);
  EXPECT_TRUE(domain.belong(point));

  const data::vector shortest = domain.shortest(data::vector(1.8, -3.5, 0.25));
  EXPECT_NEAR(shortest.x(), -0.2, 1e-12);
  EXPECT_NEAR(shortest.y(), 0.5, 1e-12);
  EXPECT_NEAR(shortest.z(), 0.25, 1e-12);
}

TEST(DomainTest, NearestDistanceDetectsPeriodicBoundaryProximity) {
  lattice::Domain domain(-1.0, 1.0, -1.0, 1.0, 0.0, 0.0);
  domain.set_periodic(0, true);

  EXPECT_DOUBLE_EQ(domain.nearest_distance(data::vector(-0.8, 0.0, 0.0), 0), 0.2);
  EXPECT_EQ(domain.nearest_distance(data::vector(-0.8, 0.0, 0.0), 1), std::numeric_limits<double>::max());
}

}  // namespace mdcraft::test

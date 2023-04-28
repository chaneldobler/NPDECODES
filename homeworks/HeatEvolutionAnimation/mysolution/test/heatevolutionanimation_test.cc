/**
 * @file heatevolutionanimation_test.cc
 * @brief NPDE homework HeatEvolutionAnimation code
 * @author Oliver Rietmann, Erick Schulz
 * @date 01.01.2020
 * @copyright Developed at ETH Zurich
 */

#include "../heatevolutionanimation.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

namespace HeatEvolutionAnimation::test {

TEST(HeatEvolutionAnimation, dummyFunction) {
  double x = 0.0;
  int n = 0;

  Eigen::Vector2d v = HeatEvolutionAnimation::dummyFunction(x, n);

  Eigen::Vector2d v_ref = {1.0, 1.0};

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, (v - v_ref).lpNorm<Eigen::Infinity>(), tol);
}

}  // namespace HeatEvolutionAnimation::test

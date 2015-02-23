#include <gtest/gtest.h>
#include <src_jf_aligner/frag_info.hpp>
#include <src_jf_aligner/pb_aligner.hpp>

namespace {
TEST(ComputeCoordsInfo, LeastSquare1) {
  std::vector<int> ul({50});
  frag_lists::frag_info frag(50, "1F");
  align_pb::mer_lists   ml(&frag);
  ml.fwd.offsets = {{10, 20}};
  ml.fwd.lis     = {0};

  auto res = compute_coords_info(ml, 50, 13, 70, &ul, true);
  EXPECT_EQ(1.0, res.stretch);
  EXPECT_EQ(-10.0, res.offset);
  EXPECT_EQ(0.0, res.avg_err);
}

TEST(ComputeCoordsInfo, LeastSquareConsecutive) {
  std::vector<int> ul({50});
  frag_lists::frag_info frag(50, "1F");
  align_pb::mer_lists   ml(&frag);
  ml.fwd.offsets = {{10, 20}, {11, 21}, {12, 22}, {13, 23}};
  ml.fwd.lis     = {0, 1, 2, 3};

  auto res = compute_coords_info(ml, 50, 13, 70, &ul, true);
  EXPECT_NEAR(1.0, res.stretch, 1e-6);
  EXPECT_NEAR(-10.0, res.offset, 1e-6);
  EXPECT_NEAR(0.0, res.avg_err, 1e-6);
}
} // empty namespace

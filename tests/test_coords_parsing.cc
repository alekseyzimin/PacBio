#include <gtest/gtest.h>
#include <src_jf_aligner/coords_parsing.hpp>

#include <gtest/gtest.h>

namespace {
TEST(CoordsParsing, CoordsInfo) {
  std::istringstream is("339 452 33 137 22 19 20 64 65 8520 137 1.11185 302.084 0.727273 116700F_68092F 13:41 6:34 15:57");
  align_pb::coords_info ci;
  is >> ci;

  EXPECT_EQ(339, ci.rs);
  EXPECT_EQ(452, ci.re);
  EXPECT_EQ(33, ci.qs);
  EXPECT_EQ(137, ci.qe);
  EXPECT_EQ(22, ci.nb_mers);
  EXPECT_EQ((unsigned int)19, ci.pb_cons);
  EXPECT_EQ((unsigned int)20, ci.sr_cons);
  EXPECT_EQ((unsigned int)64, ci.pb_cover);
  EXPECT_EQ((unsigned int)65, ci.sr_cover);
  EXPECT_EQ((size_t)8520, ci.rl);
  EXPECT_EQ((size_t)137, ci.ql);
  EXPECT_DOUBLE_EQ(1.11185, ci.stretch);
  EXPECT_DOUBLE_EQ(302.084, ci.offset);
  EXPECT_DOUBLE_EQ(0.727273, ci.avg_err);
  EXPECT_EQ("116700F_68092F", ci.qname);
  EXPECT_EQ("116700F_68092F", ci.unitigs.name());
  EXPECT_EQ(std::vector<int>({13, 6, 15}), ci.kmers_info);
  EXPECT_EQ(std::vector<int>({41, 34, 57}), ci.bases_info);
}
} // empty namespace

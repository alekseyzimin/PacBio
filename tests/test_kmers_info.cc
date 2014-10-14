#include <gtest/gtest.h>

#include <src_jf_aligner/pb_aligner.hpp>

namespace {
typedef std::vector<int> vi;
struct mock_align_pb {
  std::unique_ptr<vi> unitigs_lengths_;
  unsigned int        k_len_;
};

TEST(ComputeKmersInfo, SimpleOverlap) {
  mer_dna::k(17);
  static const unsigned int k_len = 31;
  const std::vector<int> ul({ 100, 100, 100 });
  //  aligner.unitigs_lengths_.reset(new vi({ 100, 100, 100 }));
  const std::string bad_name = "0F_1R_3F";
  const std::string good_name = "0F_1R_2F";
  const super_read_name bad_sr(bad_name);
  const super_read_name good_sr(good_name);

  vi bad_mer_info, good_mer_info;
  vi bad_base_info, good_base_info;
  align_pb::compute_kmers_info compute_bad(bad_mer_info, bad_base_info, bad_sr, k_len, mer_dna::k(), &ul);
  align_pb::compute_kmers_info compute_good(good_mer_info, good_base_info, good_sr, k_len, mer_dna::k(),  &ul);
  EXPECT_EQ(vi({0, 0, 0, 0, 0}), bad_mer_info);
  EXPECT_EQ(vi({0, 0, 0, 0, 0}), good_mer_info);
  EXPECT_EQ(vi({0, 0, 0, 0, 0}), bad_base_info);
  EXPECT_EQ(vi({0, 0, 0, 0, 0}), good_base_info);

  compute_bad.add_mer(20);
  EXPECT_EQ(vi({1, 0, 0, 0, 0}), bad_mer_info);
  EXPECT_EQ(vi({17, 0, 0, 0, 0}), bad_base_info);
  compute_bad.add_mer(71);
  EXPECT_EQ(vi({2, 1, 1, 0, 0}), bad_mer_info);
  EXPECT_EQ(vi({34, 17, 17, 0, 0}), bad_base_info);
  compute_bad.add_mer(85);
  EXPECT_EQ(vi({2, 1, 2, 0, 0}), bad_mer_info);
  EXPECT_EQ(vi({47, 30, 31, 0, 0}), bad_base_info);
  compute_bad.add_mer(142);
  EXPECT_EQ(vi({}), bad_mer_info);
  EXPECT_EQ(vi({}), bad_base_info);
  compute_bad.add_mer(170);
  EXPECT_EQ(vi({}), bad_mer_info);
  EXPECT_EQ(vi({}), bad_base_info);

  compute_good.add_mer(70);
  EXPECT_EQ(vi({1, 0, 0, 0, 0}), good_mer_info);
  EXPECT_EQ(vi({17, 16, 16, 0, 0}), good_base_info);
  compute_good.add_mer(84);
  EXPECT_EQ(vi({2, 1, 1, 0, 0}), good_mer_info);
  EXPECT_EQ(vi({31, 30, 30, 0, 0}), good_base_info);
  compute_good.add_mer(130);
  EXPECT_EQ(vi({2, 1, 2, 0, 0}), good_mer_info);
  EXPECT_EQ(vi({31, 30, 47, 6, 6}), good_base_info);
  compute_good.add_mer(150);
  EXPECT_EQ(vi({2, 1, 3, 1, 1}), good_mer_info);
  EXPECT_EQ(vi({31, 30, 64, 23, 23}), good_base_info);
  compute_good.add_mer(165);
  EXPECT_EQ(vi({2, 1, 3, 1, 2}), good_mer_info);
  EXPECT_EQ(vi({31, 30, 68, 27, 38}), good_base_info);
} // PbAligner.ComputeKmersInfo

TEST(ComputeKmersInfo, ComplexOverlap) {
  mer_dna::k(17);
  static const unsigned int k_len = 31;
  std::vector<int> unitigs_lengths({ 100, 31, 31, 40, 100 });
  const std::string name = "0F_1R_2F_3R_4F";
  const super_read_name sr(name);

  vi mer_info, base_info;
  align_pb::compute_kmers_info compute(mer_info, base_info, sr, k_len, mer_dna::k(), &unitigs_lengths);
  EXPECT_EQ(vi({0, 0, 0, 0, 0, 0, 0, 0, 0}), mer_info);
  EXPECT_EQ(vi({0, 0, 0, 0, 0, 0, 0, 0, 0}), base_info);

  compute.add_mer(70);
  EXPECT_EQ(vi({1, 0, 0, 0, 0, 0, 0, 0, 0}), mer_info);
  EXPECT_EQ(vi({17, 16, 16, 15, 15, 14, 14, 4, 4}), base_info);
  compute.add_mer(71);
  EXPECT_EQ(vi({2, 1, 1, 0, 0, 0, 0, 0, 0}), mer_info);
  EXPECT_EQ(vi({18, 17, 17, 16, 16, 15, 15, 5, 5}), base_info);
  compute.add_mer(72);
  EXPECT_EQ(vi({3, 2, 2, 1, 1, 0, 0, 0, 0}), mer_info);
  EXPECT_EQ(vi({19, 18, 18, 17, 17, 16, 16, 6, 6}), base_info);
  compute.add_mer(73);
  EXPECT_EQ(vi({4, 3, 3, 2, 2, 1, 1, 0, 0}), mer_info);
  EXPECT_EQ(vi({20, 19, 19, 18, 18, 17, 17, 7, 7}), base_info);
  compute.add_mer(74);
  EXPECT_EQ(vi({5, 4, 4, 3, 3, 2, 2, 0, 0}), mer_info);
  EXPECT_EQ(vi({21, 20, 20, 19, 19, 18, 18, 8, 8}), base_info);
  compute.add_mer(82);
  EXPECT_EQ(vi({6, 5, 5, 4, 4, 3, 3, 0, 0}), mer_info);
  EXPECT_EQ(vi({29, 28, 28, 27, 27, 26, 26, 16, 16}), base_info);
  compute.add_mer(83);
  EXPECT_EQ(vi({7, 6, 6, 5, 5, 4, 4, 1, 1}), mer_info);
  EXPECT_EQ(vi({30, 29, 29, 28, 28, 27, 27, 17, 17}), base_info);
  compute.add_mer(84);
  EXPECT_EQ(vi({8, 7, 7, 6, 6, 5, 5, 2, 2}), mer_info);
  EXPECT_EQ(vi({31, 30, 30, 29, 29, 28, 28, 18, 18}), base_info);
  compute.add_mer(85);
  EXPECT_EQ(vi({8, 7, 8, 7, 7, 6, 6, 3, 3}), mer_info);
  EXPECT_EQ(vi({31, 30, 31, 30, 30, 29, 29, 19, 19}), base_info);
  compute.add_mer(86);
  EXPECT_EQ(vi({8, 7, 8, 7, 8, 7, 7, 4, 4}), mer_info);
  EXPECT_EQ(vi({31, 30, 31, 30, 31, 30, 30, 20, 20}), base_info);
  compute.add_mer(87);
  EXPECT_EQ(vi({8, 7, 8, 7, 8, 7, 8, 5, 5}), mer_info);
  EXPECT_EQ(vi({31, 30, 31, 30, 31, 30, 31, 21, 21}), base_info);
  compute.add_mer(96);
  EXPECT_EQ(vi({8, 7, 8, 7, 8, 7, 9, 6, 6}), mer_info);
  EXPECT_EQ(vi({31, 30, 31, 30, 31, 30, 40, 30, 30}), base_info);
  compute.add_mer(97);
  EXPECT_EQ(vi({8, 7, 8, 7, 8, 7, 9, 6, 7}), mer_info);
  EXPECT_EQ(vi({31, 30, 31, 30, 31, 30, 40, 30, 31}), base_info);
  compute.add_mer(166);
  EXPECT_EQ(vi({8, 7, 8, 7, 8, 7, 9, 6, 8}), mer_info);
  EXPECT_EQ(vi({31, 30, 31, 30, 31, 30, 40, 30, 48}), base_info);
  compute.add_mer(167);
  EXPECT_EQ(vi({}), mer_info);
  EXPECT_EQ(vi({}), base_info);
} // ComputeKmersInfo.ComplexOverlap

} // namespace

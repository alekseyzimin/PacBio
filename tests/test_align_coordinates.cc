#include <fstream>
#include <string>
#include <gtest/gtest.h>
#include <src_jf_aligner/superread_parser.hpp>
#include <src_jf_aligner/pb_aligner.hpp>

namespace {
class FragsCoords : public ::testing::Test {
protected:
  virtual void SetUp() {
    std::ifstream file(sr_file);
    file.seekg(0, std::ios::end);
    super_read_approx_len = file.tellg();

    std::string line;
    {
      std::ifstream pb(pb_file);
      pb.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      while(std::getline(pb, line))
        pb_sequence += line;
    }

    {
      std::ifstream sr(sr_file);
      std::string* seq;
      while(sr.peek() != EOF) {
        std::getline(sr, line);
        if(line[0] == '>')
          seq = &sr_sequences[line.substr(1)];
        else
          *seq += line;
      }
    }
  }

  size_t                             super_read_approx_len;
  std::string                        pb_sequence;
  std::map<std::string, std::string> sr_sequences;
  static constexpr const char* pb_file           = "test_pacbio.fa";
  static constexpr const char* sr_file           = "test_super_reads.fa";
  static constexpr const int unitig_lengths_[8] = {
    0, 3000, 1043, 733, 1043, 1044, 0, 1822
  };
};

char rev_comp(const char c) {
  switch(c) {
  case 'A': case 'a': return 'T';
  case 'C': case 'c': return 'G';
  case 'G': case 'g': return 'C';
  case 'T': case 't': return 'A';
  }
  return 'N';
}
std::string rev_comp(const std::string s) {
  std::string res;
  for(auto it = s.crbegin(); it != s.crend(); ++it)
    res += rev_comp(*it);
  return res;
}

TEST_F(FragsCoords, Normal) {
  mer_dna::k(17);
  mer_pos_hash_type hash(super_read_approx_len * 2);
  frag_lists names(1);
  superread_parse(1, hash, names, sr_file);
  align_pb aligner(hash, 10, 2);
  auto res = aligner.align_sequence(pb_sequence);

  auto& frags_pos = res.second;
  EXPECT_EQ((size_t)3, frags_pos.size());

  // Check that all k-mers have equal sequence and that the LIS are in
  // increasing order
  for(auto& it : frags_pos) {
    auto&                   sr_seq = sr_sequences[std::string(it.first)];
    auto&                   ml     = it.second;
    for(auto offsets : ml.fwd.offsets)
      EXPECT_EQ(pb_sequence.substr(offsets.first - 1, mer_dna::k()),
                sr_seq.substr(offsets.second - 1, mer_dna::k()));

    align_pb::pb_sr_offsets p_offsets(0, 0);
    for(auto lis_id : ml.fwd.lis) {
      EXPECT_LT(p_offsets.first, ml.fwd.offsets[lis_id].first);
      EXPECT_LT(p_offsets.second, ml.fwd.offsets[lis_id].second);
      p_offsets = ml.fwd.offsets[lis_id];
    }

    for(auto offsets : ml.bwd.offsets)
      EXPECT_EQ(pb_sequence.substr(offsets.first - 1, mer_dna::k()),
                rev_comp(sr_seq.substr(-offsets.second - 1, mer_dna::k())));

    p_offsets = align_pb::pb_sr_offsets(0, std::numeric_limits<int>::min());
    for(auto lis_id : ml.bwd.lis) {
      EXPECT_LT(p_offsets.first, ml.bwd.offsets[lis_id].first);
      EXPECT_LT(p_offsets.second, ml.bwd.offsets[lis_id].second);
      p_offsets = ml.bwd.offsets[lis_id];
    }


    // But that the one before the first and after the last are different
    EXPECT_NE(pb_sequence.substr(ml.fwd.offsets.front().first - 2, mer_dna::k()),
              sr_seq.substr(ml.fwd.offsets.front().second - 2, mer_dna::k()));
    EXPECT_NE(pb_sequence.substr(ml.fwd.offsets.back().first, mer_dna::k()),
              sr_seq.substr(ml.fwd.offsets.back().second, mer_dna::k()));

    EXPECT_NE(pb_sequence.substr(ml.bwd.offsets.front().first - 2, mer_dna::k()),
              sr_seq.substr(-ml.bwd.offsets.front().second - 2, mer_dna::k()));
    EXPECT_NE(pb_sequence.substr(ml.bwd.offsets.back().first, mer_dna::k()),
              rev_comp(sr_seq.substr(-ml.bwd.offsets.back().second, mer_dna::k())));
  }
} // ComputeCoords.Consistency

} // namespace

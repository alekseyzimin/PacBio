#include <fstream>
#include <string>
#include <map>
#include <algorithm>
#include <random>
#include <gtest/gtest.h>
#include <src_jf_aligner/superread_parser.hpp>
#include <src_jf_aligner/pb_aligner.hpp>
#include <tests/misc.hpp>

namespace {
class FragsCoords : public ::testing::Test {
protected:
  virtual void SetUp() {
    std::string line;
    {
      std::ifstream lengths(ul_file);
      int i, l;
      while(true) {
        lengths >> i >> l;
        if(!lengths.good()) break;
        unitigs_lengths.push_back(l);
      }
    }

    {
      std::ifstream file(sr_file);
      file.seekg(0, std::ios::end);
      super_read_approx_len = file.tellg();
    }

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

  void check_mers_sequence(const align_pb::frags_pos_type& frags_pos) const;

  size_t                             super_read_approx_len;
  std::string                        pb_sequence;
  std::map<std::string, std::string> sr_sequences;
  static constexpr const char* pb_file           = "aligner_output/test_pacbio.fa";
  static constexpr const char* sr_file           = "aligner_output/test_super_reads.fa";
  static constexpr const char* ul_file           = "aligner_output/test_unitigs_lengths";
  std::vector<int> unitigs_lengths;
  static const int mer_len = 65; // mer len for k-unitigs
};

void FragsCoords::check_mers_sequence(const align_pb::frags_pos_type& frags_pos) const {
  // Check that all k-mers have equal sequence and that the LIS are in
  // increasing order
  for(auto& it : frags_pos) {
    auto        sr_seq_it = sr_sequences.find(std::string(it.first));
    ASSERT_NE(sr_sequences.end(), sr_seq_it);
    const auto& sr_seq    = sr_seq_it->second;
    const auto& ml        = it.second;
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
                misc::rev_comp(sr_seq.substr(-offsets.second - 1, mer_dna::k())));

    p_offsets = align_pb::pb_sr_offsets(0, std::numeric_limits<int>::min());
    for(auto lis_id : ml.bwd.lis) {
      EXPECT_LT(p_offsets.first, ml.bwd.offsets[lis_id].first);
      EXPECT_LT(p_offsets.second, ml.bwd.offsets[lis_id].second);
      p_offsets = ml.bwd.offsets[lis_id];
    }


    // But that the k-mer before the first and after the last are different
    EXPECT_NE(pb_sequence.substr(ml.fwd.offsets.front().first - 2, mer_dna::k()),
              sr_seq.substr(ml.fwd.offsets.front().second - 2, mer_dna::k()));
    EXPECT_NE(pb_sequence.substr(ml.fwd.offsets.back().first, mer_dna::k()),
              sr_seq.substr(ml.fwd.offsets.back().second, mer_dna::k()));

    EXPECT_NE(pb_sequence.substr(ml.bwd.offsets.front().first - 2, mer_dna::k()),
              sr_seq.substr(-ml.bwd.offsets.front().second - 2, mer_dna::k()));
    EXPECT_NE(pb_sequence.substr(ml.bwd.offsets.back().first, mer_dna::k()),
              misc::rev_comp(sr_seq.substr(-ml.bwd.offsets.back().second, mer_dna::k())));
  }
}

// void check_info(const align_pb::coords_info e, const align_pb::coords_info a) {
//   EXPECT_EQ(e.rs, a.rs);
//   EXPECT_EQ(e.re, a.re);
//   EXPECT_EQ(e.qs, a.qs);
//   EXPECT_EQ(e.qe, a.qe);
//   EXPECT_EQ(e.nb_mers, a.nb_mers);
//   EXPECT_EQ(e.pb_cons, a.pb_cons);
//   EXPECT_EQ(e.sr_cons, a.sr_cons);
//   EXPECT_EQ(e.pb_cover, a.pb_cover);
//   EXPECT_EQ(e.sr_cover, a.sr_cover);
//   EXPECT_EQ(e.rl, a.rl);
//   EXPECT_EQ(e.ql, a.ql);
//   EXPECT_EQ(e.rn, a.rn);
//   EXPECT_EQ(e.qname, a.unitigs.name());
//   EXPECT_NEAR(e.stretch, a.stretch, 1e-4);
//   EXPECT_NEAR(e.offset, a.offset, 1e-1);
//   EXPECT_NEAR(e.avg_err, a.avg_err, 1e-3);
// }

// void check_coords(std::map<std::string, align_pb::coords_info> res, const align_pb::coords_info_type& coords) {
//   for(const auto& info : coords) {
//     auto it = res.find(info.unitigs.name());
//     ASSERT_NE(res.end(), it);
//     check_info(it->second, info);
//   }
// }

// TEST_F(FragsCoords, NormalConsistency) {
//   mer_dna::k(17);
//   mer_pos_hash_type hash(super_read_approx_len * 2);
//   frag_lists names(1);
//   superread_parse(1, hash, names, sr_file);
//   align_pb aligner(hash, 10, 2, 100, 1);
//   auto res = aligner.align_sequence(pb_sequence);

//   auto& frags_pos = res.second;
//   EXPECT_EQ((size_t)3, frags_pos.size());
//   check_mers_sequence(frags_pos);

//   std::map<std::string, align_pb::coords_info> align_res;
//   align_res["1R_3F"] = align_pb::coords_info(302, 876, 3051, 3597, 90, 79, 79, 263, 262, 1287, 3668, false,
//                                              1.05042, -2899.38, 1.72222, "1R_3F");
//   align_res["5F_4R_2F"] = align_pb::coords_info(302, 1149, 1420, 613, 107, 93, 93, 328, 327, 1287, 3000, false,
//                                                 -1.04824, 1793.59, 1.4486, "5F_4R_2F");
//   align_res["7R_2F"] = align_pb::coords_info(302, 1149, 1280, 473, 107, 93, 93, 328, 327, 1287, 2800, false,
//                                              -1.04824, 1646.84, 1.4486, "7R_2F");
//   auto& coords = res.first;
//   EXPECT_EQ((size_t)3, coords.size());
//   check_coords(align_res, coords);
// } // ComputeCoords.Consistency

// TEST_F(FragsCoords, ForwardConsistency) {
//   mer_dna::k(17);
//   mer_pos_hash_type hash(super_read_approx_len * 2);
//   frag_lists names(1);
//   superread_parse(1, hash, names, sr_file);
//   align_pb aligner(hash, 10, 2, 20, 1, true);
//   aligner.unitigs_lengths(&unitigs_lengths, mer_len);
//   auto res = aligner.align_sequence(pb_sequence);
//   auto res_max = aligner.align_sequence(pb_sequence);

//   auto& frags_pos = res.second;
//   EXPECT_EQ((size_t)3, frags_pos.size());
//   check_mers_sequence(frags_pos);

//   auto& coords = res.first;
//   EXPECT_EQ((size_t)3, coords.size());
//   // Check that the alternate sums add up to the right values
//   for(const auto& info : coords) {
//     int total_mers  = 0;
//     int total_bases = 0;
//     ASSERT_EQ(info.kmers_info.size(), info.bases_info.size());
//     for(size_t i = 0; i < info.kmers_info.size(); ++i) {
//       total_mers  += i & 0x1 ? -info.kmers_info[i] : info.kmers_info[i];
//       total_bases += i & 0x1 ? -info.bases_info[i] : info.bases_info[i];
//     }
//     EXPECT_EQ(info.nb_mers, total_mers);
//     EXPECT_EQ((int)info.sr_cover, total_bases);
//   }

//   std::map<std::string, align_pb::coords_info> align_res;
//   align_res["1R_3F"] = align_pb::coords_info(302, 876, 3051, 3597, 90, 79, 79, 263, 262, 1287, 3668, false,
//                                              1.05042, -2899.38, 1.72222, "1R_3F");
//   align_res["2R_4F_5R"] = align_pb::coords_info(302, 1149, 1581, 2388, 107, 93, 93, 328, 327, 1287, 3000, true,
//                                                 1.04824, -1351.17, 1.4486, "2R_4F_5R");
//   align_res["2R_7F"] = align_pb::coords_info(302, 1149, 1521, 2328, 107, 93, 93, 328, 327, 1287, 2800, true,
//                                              1.04824, -1288.28, 1.4486, "2R_7F");

//   check_coords(align_res, coords);

//   auto& coords_max = res_max.first;
//   check_coords(align_res, coords_max);

//   // Same but with filtering
//   {
//     align_pb aligner_mer_filter(hash, 10, 2, 20, 1, true, false, 0, 0.09);
//     auto res_mer_filter = aligner_mer_filter.align_sequence(pb_sequence);
//     auto res_max_mer_filter = aligner_mer_filter.align_sequence_max(pb_sequence);
//     auto& coords_mer_filter = res_mer_filter.first;
//     auto& coords_max_mer_filter = res_max_mer_filter.first;
//     ASSERT_EQ((size_t)1, coords_mer_filter.size());
//     ASSERT_EQ((size_t)1, coords_max_mer_filter.size());
//     EXPECT_EQ("1R_3F", coords_mer_filter[0].qname);
//     EXPECT_EQ("1R_3F", coords_max_mer_filter[0].qname);
//     check_coords(align_res, coords_mer_filter);
//     check_coords(align_res, coords_max_mer_filter);
//   }

//   {
//     align_pb aligner_mer_filter(hash, 10, 2, 20, 1, true, false, 0, 0.0, 0.27);
//     auto res_mer_filter = aligner_mer_filter.align_sequence(pb_sequence);
//     auto res_max_mer_filter = aligner_mer_filter.align_sequence_max(pb_sequence);
//     auto& coords_mer_filter = res_mer_filter.first;
//     auto& coords_max_mer_filter = res_max_mer_filter.first;
//     ASSERT_EQ((size_t)1, coords_mer_filter.size());
//     ASSERT_EQ((size_t)1, coords_max_mer_filter.size());
//     EXPECT_EQ("1R_3F", coords_mer_filter[0].qname);
//     EXPECT_EQ("1R_3F", coords_max_mer_filter[0].qname);
//     check_coords(align_res, coords_mer_filter);
//     check_coords(align_res, coords_max_mer_filter);
//   }
// }

TEST(DiscardLIS, Discard) {
  align_pb::off_lis                    ol;
  std::vector<unsigned int>            all_indices;
  std::vector<align_pb::pb_sr_offsets> all_offsets;
  const int                            max_size = 100;
  auto                                 rng      = std::default_random_engine();
  for(int i = 0; i < max_size; ++i) {
    all_offsets.push_back(std::make_pair(i, 5 * i));
    all_indices.push_back(i);
  }

  // Random lis of increasing size
  for(int i = 0; i <= max_size; ++i) {
    SCOPED_TRACE(::testing::Message() << "i:" << i);
    ol.offsets = all_offsets;
    ol.lis     = all_indices;
    std::shuffle(ol.lis.begin(), ol.lis.end(), rng);
    ol.lis.resize(i);
    std::sort(ol.lis.begin(), ol.lis.end());
    ol.discard_LIS();

    ASSERT_EQ((size_t)(max_size - i), ol.offsets.size());
    int          j   = 0, k = 0;
    for(auto it : ol.lis) {
      for( ; j < (int)(it - k); ++j) {
        SCOPED_TRACE(::testing::Message() << "j:" << j);
        EXPECT_EQ(j + k, ol.offsets[j].first);
      }
      ++k;
    }
  }
}

} // namespace

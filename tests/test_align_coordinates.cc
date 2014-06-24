#include <gtest/gtest.h>

namespace {
class ComputeCoordsConsistency : public ::testing::Test {
public:
  ComputeCoordsConsistency();

protected:
  virtual void SetUp() {
    std::ifstream file(sr_file);
    file.seekg(0, std::ios::end);
    super_read_approx_len = file.tellg();
  }

  std::vector<coords_info> compute_coordinates(bool forward) {
    file_vector files;
    files.push_back(pb_path);
    stream_manager streams(files.cbegin(), files.cend());
    align_pb aligner(1 /*threads*/, hash, streams, 10 /*stretch_constant*/, 2/*stretch_factor*/,
                     forward, false /*compress*/, false /*duplicate*/,
                     0.0 /*matching_mers*/, 0.0/*matching_bases*/);
    std::vector<int> unitig_lengths.insert(unitig_lengths_);
    aligner.unitigs_lengths(&unitigs_lengths, 70 /* k for creating k-unitigs */)
      
  }

  size_t super_read_approx_len;
  static constexpr const char* pb_file = "test_pacbio.fa";
  static constexpr const char* sr_file = "test_super_reads.fa";
  static constexprt const int unitig_lengths_[8] = {
    0, 3000, 1043, 733, 1043, 1044, 0, 1822
  };
};

};

TEST_F(ComputeCoordsConsistency, Normal) {
  mer_dna::k(17);
  mer_pos_hash_type hash(super_read_approx_len * 2);
  frag_lists names(1);
  superread_parse(1, hash, names, sr_file);
  align_pb aligner(1 /*threads*/, hash, 
} // ComputeCoords.Consistency

} // namespace

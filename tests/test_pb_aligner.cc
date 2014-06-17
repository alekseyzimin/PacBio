#include <gtest/gtest.h>
#include <src_jf_aligner/jf_aligner.hpp>
#include <src_jf_aligner/pb_aligner.hpp>
#include <src_jf_aligner/superread_parser.hpp>

namespace {
char rev_char(char c) {
  switch(c) {
  case 'A': return 'T';
  case 'C': return 'G';
  case 'G': return 'C';
  case 'T': return 'A';
  }
  return 'N';
}
std::string rev_comp(std::string s) {
  for(size_t i = 0, j = s.size() - 1; i <= j; ++i, --j) {
    char tmp = rev_char(s[i]);
    s[i] = rev_char(s[j]);
    s[j] = tmp;
  }
  return s;
}

// Create a fake sequence of 1000 bp. Then a PacBio read of length 500
// starting at base 300 of the sequence. The PacBio reads has an error
// a base 100 and an insertion at base 400.
//
// A bunch of super-reads are created, all of length 100. R1 and R2
// spans the left boundary of the PacBio read, R3 and R4 are over the
// sequencing error at base 100, R5 and R6 are over the insertion of
// length 4 at base 400, R7 and R8 span the right boundary of the
// PacBio read, finally R9 and R10 are in the error free area (between
// base 101 and 399 of the PacBio read).
//
// An odd numbered read is forward, and even numbered read is reversed.
//
// This is the picture. X representing the sequencing error on I the
// insertion point.
//
//    Full sequence
// 0     300          400                         700                     1000
// -----------------------------------------------------------------------
// PacBio --------------X---------------------------I------------
//     |----> R1       |----> R3  |----> R9      |----> R5    |----> R7
//      <----| R2    <----| R4      <----| R10    <----| R6    <----| R8

std::string generate_sequences(const char* pacbio_reads, const char* superreads) {
  std::string sequence(1000, 'A');
  static const char base[4] = {'A', 'C', 'G', 'T' };
  for(size_t i = 0; i < sequence.size(); ++i)
    sequence[i] = base[random() % 4];

  std::string pacbio_sequence = sequence.substr(300, 500);
  char error = pacbio_sequence[100];
  do {
    pacbio_sequence[100] = base[random() % 4];
  } while(error == pacbio_sequence[100]);
  pacbio_sequence[399] = 'C';
  pacbio_sequence[400] = 'G';
  pacbio_sequence.insert(400, "ACGT");
  {
    std::ofstream pacbio(pacbio_reads);
    if(!pacbio.good())
      throw std::runtime_error("Can't open pacbio temp file");
    pacbio << ">pacbio\n" << pacbio_sequence << "\n";
  }

  {
    std::ofstream super(superreads);
    if(!super.good())
      throw std::runtime_error("Can't open superreads temp file");
    super << ">R1\n" << sequence.substr(250, 100) << "\n"
          << ">R2\n" << rev_comp(sequence.substr(275, 100)) << "\n"
          << ">R3\n" << sequence.substr(390, 100) << "\n"
          << ">R4\n" << rev_comp(sequence.substr(380, 100)) << "\n"
          << ">R5\n" << sequence.substr(660, 100) << "\n"
          << ">R6\n" << rev_comp(sequence.substr(670, 100)) << "\n"
          << ">R7\n" << sequence.substr(770, 100) << "\n"
          << ">R8\n" << rev_comp(sequence.substr(780, 100)) << "\n"
          << ">R9\n" << sequence.substr(500, 100) << "\n"
          << ">R10\n" << rev_comp(sequence.substr(550, 100)) << "\n";
  }

  return pacbio_sequence;
}

struct remove_file {
  const char* path;
  bool do_unlink;
  remove_file(const char* p, bool unlink = true) : path(p), do_unlink(unlink) { }
  ~remove_file() { if(do_unlink) unlink(path); }
};

TEST(PbAligner, FakeSequences) {
  mer_dna::k(15);
  remove_file pb_file(".pacbio.fa", false);
  remove_file sr_file(".superreads.fa", false);
  std::string pacbio_sequence = generate_sequences(pb_file.path, sr_file.path);

  mer_pos_hash_type hash(2048);
  frag_lists names(1);
  superread_parse(1, hash, names, sr_file.path);
  EXPECT_EQ((size_t)1, names.size());

  parse_sequence parser(pacbio_sequence);
  align_pb::frags_pos_type frags_pos;
  align_pb::process_read(hash, parser, frags_pos, 10, 2);

  EXPECT_EQ((size_t)10, frags_pos.size());
  for(auto it = frags_pos.cbegin(); it != frags_pos.cend(); ++it) {
    const align_pb::mer_lists& ml = it->second;
    int read_id = std::atoi(it->first + 1);
    ASSERT_TRUE(read_id >= 1 && read_id <= 10); // Read id is valid
    EXPECT_TRUE(std::is_sorted(ml.fwd.offsets.cbegin(), ml.fwd.offsets.cend())); // mers offsets must be sorted
    EXPECT_TRUE(std::is_sorted(ml.bwd.offsets.cbegin(), ml.bwd.offsets.cend())); // mers offsets must be sorted
    EXPECT_EQ(ml.fwd.offsets.size(), ml.fwd.lis.size()); // lis has same size
    {
      auto iit = ml.fwd.lis.cbegin();
      for(auto it = ml.fwd.offsets.cbegin(); it != ml.fwd.offsets.cend(); ++it, ++iit) {
        EXPECT_EQ(it->second, ml.fwd.offsets[*iit].second); // and is equivalent
        EXPECT_EQ(it->first, ml.fwd.offsets[*iit].first); // and is equivalent
      }
    }

    EXPECT_EQ(ml.bwd.offsets.size(), ml.bwd.lis.size()); // lis has same size
    {
      auto iit = ml.bwd.lis.cbegin();
      for(auto it = ml.bwd.offsets.cbegin(); it != ml.bwd.offsets.cend(); ++it, ++iit) {
        EXPECT_EQ(it->second, ml.bwd.offsets[*iit].second); // and is equivalent
        EXPECT_EQ(it->first, ml.bwd.offsets[*iit].first); // and is equivalent
      }
    }
    SCOPED_TRACE(::testing::Message() << "Read:" << read_id);

    switch(read_id) {
    case 1:
      ASSERT_EQ((size_t)36, ml.fwd.offsets.size());
      EXPECT_EQ(51, ml.fwd.offsets.front().second);
      EXPECT_EQ(86, ml.fwd.offsets.back().second);
      break;

    case 2:
      ASSERT_EQ((size_t)61, ml.bwd.offsets.size());
      EXPECT_EQ(-61, ml.bwd.offsets.front().second);
      EXPECT_EQ(-1, ml.bwd.offsets.back().second);
      break;

    case 3:
      ASSERT_EQ((size_t)75, ml.fwd.offsets.size());
      EXPECT_EQ(12, ml.fwd.offsets.front().second);
      EXPECT_EQ(86, ml.fwd.offsets.back().second);
      break;

    case 4:
      ASSERT_EQ((size_t)71, ml.bwd.offsets.size());
      EXPECT_EQ(-86, ml.bwd.offsets.front().second);
      EXPECT_EQ(-1, ml.bwd.offsets.back().second);
      break;

    case 5:
      ASSERT_EQ((size_t)70, ml.fwd.offsets.size());
      EXPECT_EQ(1, ml.fwd.offsets.front().second);
      EXPECT_EQ(86, ml.fwd.offsets.back().second);
      break;

    case 6:
      ASSERT_EQ((size_t)70, ml.bwd.offsets.size());
      EXPECT_EQ(-86, ml.bwd.offsets.front().second);
      EXPECT_EQ(-1, ml.bwd.offsets.back().second);
      break;

    case 7:
      ASSERT_EQ((size_t)16, ml.fwd.offsets.size());
      EXPECT_EQ(1, ml.fwd.offsets.front().second);
      EXPECT_EQ(16, ml.fwd.offsets.back().second);
      break;

    case 8:
      ASSERT_EQ((size_t)6, ml.bwd.offsets.size());
      EXPECT_EQ(-86, ml.bwd.offsets.front().second);
      EXPECT_EQ(-81, ml.bwd.offsets.back().second);
      break;

    case 9:
      ASSERT_EQ((size_t)86, ml.fwd.offsets.size());
      EXPECT_EQ(1, ml.fwd.offsets.front().second);
      EXPECT_EQ(86, ml.fwd.offsets.back().second);
      break;

    case 10:
      ASSERT_EQ((size_t)86, ml.bwd.offsets.size());
      EXPECT_EQ(-86, ml.bwd.offsets.front().second);
      EXPECT_EQ(-1, ml.bwd.offsets.back().second);
      break;
    }
  }
}

TEST(PbAligner, ReverseSRName) {
  EXPECT_EQ("", align_pb::reverse_super_read_name(""));
  EXPECT_EQ("toto", align_pb::reverse_super_read_name("toto"));
  EXPECT_EQ("1F_2R_3_4R", align_pb::reverse_super_read_name("1F_2R_3_4R"));

  EXPECT_EQ("asdfR", align_pb::reverse_super_read_name("asdfF"));
  EXPECT_EQ("1234F", align_pb::reverse_super_read_name("1234R"));
  EXPECT_EQ("asdf;lkjqweropuF_1R", align_pb::reverse_super_read_name("1F_asdf;lkjqweropuR"));
  EXPECT_EQ("4R_3F_2R_1F", align_pb::reverse_super_read_name("1R_2F_3R_4F"));
} // PbAligner.ReverseSRName

typedef std::vector<int> vi;
struct mock_align_pb {
  std::unique_ptr<vi> unitigs_lengths_;
  unsigned int        k_len_;
};

TEST(ComputeKmersInfo, SimpleOverlap) {
  mer_dna::k(17);
  mock_align_pb aligner;
  aligner.k_len_ = 31;
  aligner.unitigs_lengths_.reset(new vi({ 100, 100, 100 }));
  const std::string bad_name = "0F_1R_3F";
  const std::string good_name = "0F_1R_2F";

  vi bad_mer_info, good_mer_info;
  vi bad_base_info, good_base_info;
  align_pb::compute_kmers_info<mock_align_pb> compute_bad(bad_mer_info, bad_base_info, bad_name, aligner);
  align_pb::compute_kmers_info<mock_align_pb> compute_good(good_mer_info, good_base_info, good_name, aligner);
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
  mock_align_pb aligner;
  aligner.k_len_ = 31;
  aligner.unitigs_lengths_.reset(new vi({ 100, 31, 31, 40, 100 }));
  const std::string name = "0F_1R_2F_3R_4F";

  vi mer_info, base_info;
  align_pb::compute_kmers_info<mock_align_pb> compute(mer_info, base_info, name, aligner);
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

} // namespace {

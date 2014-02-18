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
  remove_file pb_file("/tmp/pacbio.fa", false);
  remove_file sr_file("/tmp/superreads.fa", false);
  std::string pacbio_sequence = generate_sequences(pb_file.path, sr_file.path);

  mer_pos_hash_type hash(1024);
  name_lists names(1);
  superread_parse(1, hash, names, sr_file.path);
  EXPECT_EQ((size_t)1, names.size());

  parse_sequence parser(pacbio_sequence);
  align_pb::frags_pos_type frags_pos;
  align_pb::process_read(hash, parser, frags_pos);

  EXPECT_EQ((size_t)10, frags_pos.size());
  for(auto it = frags_pos.cbegin(); it != frags_pos.cend(); ++it) {
    const align_pb::mer_lists& ml = it->second;
    int read_id = std::atoi(it->first + 1);
    ASSERT_TRUE(read_id >= 1 && read_id <= 10); // Read id is valid
    EXPECT_TRUE(std::is_sorted(ml.offsets.cbegin(), ml.offsets.cend())); // mers offsets must be sorted
    EXPECT_EQ(ml.offsets.size(), ml.lis.size()); // lis has same size
    auto iit = ml.lis.cbegin();
    for(auto it = ml.offsets.cbegin(); it != ml.offsets.cend(); ++it, ++iit) {
      EXPECT_EQ(it->sr_offset, ml.offsets[*iit].sr_offset); // and is equivalent
      EXPECT_EQ(it->pb_offset, ml.offsets[*iit].pb_offset); // and is equivalent
    }
    SCOPED_TRACE(::testing::Message() << "Read:" << read_id);

    switch(read_id) {
    case 1:
      EXPECT_EQ((size_t)36, ml.offsets.size());
      EXPECT_EQ(50, ml.offsets.front().sr_offset);
      EXPECT_EQ(85, ml.offsets.back().sr_offset);
      break;

    case 2:
      EXPECT_EQ((size_t)61, ml.offsets.size());
      EXPECT_EQ(-60, ml.offsets.front().sr_offset);
      EXPECT_EQ(0, ml.offsets.back().sr_offset);
      break;

    case 3:
      EXPECT_EQ((size_t)75, ml.offsets.size());
      EXPECT_EQ(11, ml.offsets.front().sr_offset);
      EXPECT_EQ(85, ml.offsets.back().sr_offset);
      break;

    case 4:
      EXPECT_EQ((size_t)71, ml.offsets.size());
      EXPECT_EQ(-85, ml.offsets.front().sr_offset);
      EXPECT_EQ(0, ml.offsets.back().sr_offset);
      break;

    case 5:
      EXPECT_EQ((size_t)70, ml.offsets.size());
      EXPECT_EQ(0, ml.offsets.front().sr_offset);
      EXPECT_EQ(85, ml.offsets.back().sr_offset);
      break;

    case 6:
      EXPECT_EQ((size_t)70, ml.offsets.size());
      EXPECT_EQ(-85, ml.offsets.front().sr_offset);
      EXPECT_EQ(0, ml.offsets.back().sr_offset);
      break;

    case 7:
      EXPECT_EQ((size_t)16, ml.offsets.size());
      EXPECT_EQ(0, ml.offsets.front().sr_offset);
      EXPECT_EQ(15, ml.offsets.back().sr_offset);
      break;

    case 8:
      EXPECT_EQ((size_t)6, ml.offsets.size());
      EXPECT_EQ(-85, ml.offsets.front().sr_offset);
      EXPECT_EQ(-80, ml.offsets.back().sr_offset);
      break;

    case 9:
      EXPECT_EQ((size_t)86, ml.offsets.size());
      EXPECT_EQ(0, ml.offsets.front().sr_offset);
      EXPECT_EQ(85, ml.offsets.back().sr_offset);
      break;

    case 10:
      EXPECT_EQ((size_t)86, ml.offsets.size());
      EXPECT_EQ(-85, ml.offsets.front().sr_offset);
      EXPECT_EQ(0, ml.offsets.back().sr_offset);
      break;
    }
  }
}

} // namespace {
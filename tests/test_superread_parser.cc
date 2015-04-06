#include <unistd.h>
#include <cstdlib>
#include <string>

#include <gtest/gtest.h>
#include <src_jf_aligner/superread_parser.hpp>
#include <tests/misc.hpp>

namespace {
using misc::remove_file;

void create_sequence(const char* path) {
  static const char bases[4] = { 'A', 'C', 'G','T' };
  std::ofstream f(path);
  if(!f.good())
    throw std::runtime_error("Failed to open sequence file");
  f << ">superread\n";
  for(int i = 0; i < 100; ++i)
    f << bases[random() % 4];
  f << "\n";
}

TEST(SuperReadParser, OneRead) {
  remove_file file;
  create_sequence(file.path);

  mer_dna::k(17);
  sequence_psa psa;
  psa.append_fasta(file.path);
  psa.compute_psa(10, mer_dna::k(), 1);
  EXPECT_TRUE(psa.check_psa());

  // Check every k-mer in the sequence
  std::ifstream f(file.path);
  if(!f.good())
    throw std::runtime_error("Failed to open sequence file");
  std::string line;
  std::getline(f, line); // Skip header
  std::getline(f, line);

  for(size_t i = 0; i < line.size() - mer_dna::k() + 1; ++i) {
    mer_dna m(line.substr(i, mer_dna::k()));
    mer_dna rm(m.get_reverse_complement());
    bool is_canonical = m < rm;

    const auto list = psa.equal_range(m);
    ASSERT_TRUE(list.first != list.second);
    EXPECT_EQ(1, std::distance(list.first, list.second));
    //    const mer_pos_hash_type::position_type& pos = *list->cbegin();
    const auto it = list.first;
    EXPECT_EQ("superread", it->frag->fwd.name);
    EXPECT_EQ((int)(i + 1) * (is_canonical ? 1 : -1), it->offset);
  }
}

std::string create_sequences(const char* path, size_t size, int nb_reads, size_t read_size) {
  static const char bases[4] = { 'A', 'C', 'G','T' };
  std::string sequence(size, 'A');
  for(size_t i = 0; i < sequence.size(); ++i)
    sequence[i] = bases[random() % 4];

  std::ofstream f(path);
  if(!f.good())
    throw std::runtime_error("Failed to open sequence file");
  for(int i = 0; i < nb_reads; ++i) {
    f << ">" << i << "\n"
      << sequence.substr(i * (size - read_size) / (nb_reads - 1), read_size) << "\n";
  }
  return sequence;
}

template<typename T>
T ceil_div(T x, T y) {
  return x / y + (x % y != 0);
}

TEST(SuperReadParser, ManyReads) {
  remove_file file;
  static const int    nb_reads  = 21;
  static const int    seq_len   = 1000;
  static const int    read_len  = 100;
  static const int    delta     = (seq_len - read_len) / (nb_reads - 1);
  static const size_t total_len = nb_reads * read_len;
  const std::string   sequence  = create_sequences(file.path, seq_len, nb_reads, read_len);

  EXPECT_EQ((size_t)seq_len, sequence.size());

  mer_dna::k(17);
  static const unsigned int min_size = 10;

  sequence_psa psa;
  psa.append_fasta(file.path);
  EXPECT_EQ(total_len, psa.sequence_size());
  psa.compute_psa(min_size, mer_dna::k());
  EXPECT_EQ(total_len - (min_size - 1) * nb_reads, psa.nb_mers());
  EXPECT_TRUE(psa.check_psa());
  EXPECT_EQ((size_t)nb_reads, psa.m_headers.size());

  for(size_t i = 0; i < sequence.size() - mer_dna::k() + 1; ++i) {
    mer_dna m(sequence.substr(i, mer_dna::k()));
    mer_dna rm(m.get_reverse_complement());
    const bool is_canonical = m < rm;
    SCOPED_TRACE(::testing::Message() << "i:" << i << " m:" << m << " canonical:" << is_canonical << " delta:" << delta);

    auto list = psa.equal_range(m, rm);
    EXPECT_TRUE(list.first != list.second);
    int count = 0;
    for(auto it = list.first; it != list.second; ++it, ++count) {
      EXPECT_EQ(read_len, (int)it->frag->len);
      int read_id = std::atoi(it->frag->fwd.name.c_str());
      // Is id valid?
      EXPECT_TRUE(read_id >= 0 && read_id < nb_reads);
      // Is the read covering position i?
      EXPECT_TRUE((size_t)(read_id * delta) <= i && (size_t)(read_id * delta + read_len) > i);
      // Is offset valid
      EXPECT_EQ((int)i + 1 - read_id * delta, is_canonical ? it->offset : -it->offset);
    }
    // Is number of reads covering position i correct?
    if(i <= read_len - mer_dna::k())
      EXPECT_EQ((int)i / delta + 1, count);
    else if(i >= seq_len - read_len)
      EXPECT_EQ((seq_len - read_len) / delta - (int)(i - read_len + mer_dna::k() - 1) / delta, count);
    else
      EXPECT_EQ((int)i / delta - (int)(i - read_len + mer_dna::k() - 1) / delta, count);
  }
}
}

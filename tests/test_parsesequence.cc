#include <string>
#include <random>

#include <gtest/gtest.h>

#include <src_jf_aligner/jf_aligner.hpp>

namespace {
std::string create_sequence(int len) {
  std::default_random_engine         gen;
  std::uniform_int_distribution<int> gen4(0, 3);

  static const char bases[4] = { 'A', 'C', 'G','T' };
  std::string res;
  for(int i = 0; i < len; ++i)
    res += bases[gen4(gen) % 4];
  return res;
}

std::string create_sequence_with_N(int len, double proba_N) {
  std::default_random_engine             gen;
  std::uniform_int_distribution<int>     gen4(0, 3);
  std::uniform_real_distribution<double> genf;

  static const char bases[4] = { 'A', 'C', 'G','T' };
  std::string res;
  for(int i = 0; i < len; ++i) {
    res += genf(gen) < proba_N ? 'N' : bases[gen4(gen) % 4];
  }
  return res;
}


TEST(ParseSequence, Uncompress) {
  static const int len = 1000;
  mer_dna::k(21);

  const std::string seq = create_sequence(len);
  parse_sequence parser(seq);
  int i = 0;
  mer_dna mer;
  for( ; parser.next(); ++i) {
    ASSERT_LE(i, (int)(len - mer_dna::k()));
    mer = seq.substr(i, mer_dna::k());
    EXPECT_TRUE(mer == parser.mer<0>().m || mer == parser.mer<0>().rm);
  }
  EXPECT_EQ((int)(len - mer_dna::k() + 1), i);
}

std::string compress(const std::string& seq) {
  std::string res;
  auto it    = seq.begin();
  char prev  = *it;
  res       += prev;
  for(++it; it != seq.end(); prev = *it, ++it) {
    if(prev != *it)
      res += *it;
  }
  return res;
}


TEST(ParseSequence, Compress) {
  static const int len = 1000;
  mer_dna::k(21);

  const std::string seq = create_sequence(len);
  const std::string cseq = compress(seq);
  const int clen = cseq.length();
  EXPECT_LE(clen, (int)seq.size());
  for(int i = 1; i < clen; ++i)
    EXPECT_NE(cseq[i - 1], cseq[i]);

  parse_sequence parser(seq, true);
  parse_sequence cparser(cseq);

  while(parser.next()) {
    ASSERT_TRUE(cparser.next());
    ASSERT_EQ(cparser.mer<0>().m, parser.mer<0>().m);
    EXPECT_EQ(cparser.mer<0>().rm, parser.mer<0>().rm);
    EXPECT_LT(cparser.seq_offset, parser.seq_offset);
  }
  ASSERT_FALSE(cparser.next());
}

typedef mer_dna long_mer_type;
typedef jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1> short_mer_type;
typedef parse_sequence_ns::parser_base<long_mer_type, short_mer_type> parse_sequence2;

TEST(ParseSequence, MultipleMers) {
  static const int len = 1000;
  long_mer_type::k(21);
  short_mer_type::k(15);

  const std::string seq = create_sequence_with_N(len, 0.05);
  parse_sequence2   parser(seq);
  long_mer_type     ml;
  short_mer_type    ms;

  for( ; parser.next(); ) {
    const int i = parser.seq_offset - 1;
    const int loffset = i - long_mer_type::k() + 1;
    const int soffset = i - short_mer_type::k() + 1;

    ASSERT_LE(i, (int)(len - short_mer_type::k()));
    if(loffset >= 0)
      EXPECT_EQ(parser.mer<0>().valid,
                std::string::npos == seq.substr(loffset, long_mer_type::k()).find('N'));
    if(soffset >= 0)
      EXPECT_EQ(parser.mer<1>().valid,
                std::string::npos == seq.substr(soffset, short_mer_type::k()).find('N'));
    if(parser.mer<0>().valid) {
      ml = seq.substr(loffset, long_mer_type::k());
      EXPECT_TRUE(ml == parser.mer<0>().m || ml == parser.mer<0>().rm);
    }
    if(parser.mer<1>().valid) {
      ms = seq.substr(soffset, short_mer_type::k());
      EXPECT_TRUE(ms == parser.mer<1>().m || ms == parser.mer<1>().rm);
    }
  }
}
} // empty namespace

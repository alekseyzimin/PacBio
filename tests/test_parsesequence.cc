#include <string>

#include <gtest/gtest.h>

#include <src_jf_aligner/jf_aligner.hpp>

namespace {
std::string create_sequence(int len) {
  static const char bases[4] = { 'A', 'C', 'G','T' };
  std::string res;
  for(int i = 0; i < len; ++i)
    res += bases[random() % 4];
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
    EXPECT_TRUE(mer == parser.m || mer == parser.rm);
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
    ASSERT_EQ(cparser.m, parser.m);
    EXPECT_EQ(cparser.rm, parser.rm);
    EXPECT_LT(cparser.offset, parser.offset);
  }
  ASSERT_FALSE(cparser.next());
}
} // empty namespace

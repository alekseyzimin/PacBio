#include <gtest/gtest.h>
#include <src_jf_aligner/superread_parser.hpp>
#include <jellyfish/mer_dna.hpp>

namespace {
using jellyfish::mer_dna;
TEST(JFMer, str_to_mer) {
  static const size_t iterations = 10;
  mer_dna::k(17);
  mer_dna m;

  for(size_t i = 0; i < iterations; ++i) {
    m.randomize();
    for(unsigned off = 0; off < mer_dna::k() - 1; ++off) {
      for(unsigned len = 1; len <= mer_dna::k() - off; ++len) {
        SCOPED_TRACE(::testing::Message() << "i:" << i << " off:" << off << " len:" << len << " m:" << m);
        uint64_t mstr = mer_sa_imp::str_to_mer(m.to_str().c_str() + off, len);
        uint64_t mdna = mer_sa_imp::str_to_mer(mer_dna_ptr<mer_dna>(m) + off, len);
        EXPECT_EQ(mstr, mdna);
      }
    }
  }
}
} // empty namespace

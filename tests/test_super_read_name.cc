#include <gtest/gtest.h>
#include <src_jf_aligner/super_read_name.hpp>
#include <tests/misc.hpp>

namespace {
TEST(SuperReadName, Parse) {
  super_read_name n1("");
  EXPECT_EQ((size_t)0, n1.nb_unitigs());
  EXPECT_EQ("", n1.get_reverse().name());
  {
    const auto u = n1[0];
    EXPECT_EQ(super_read_name::invalid_id, u.id());
    EXPECT_EQ(u.id(), n1.unitig_id(0));
  }

  super_read_name n2("1234F");
  EXPECT_EQ((size_t)1, n2.nb_unitigs());
  EXPECT_EQ("1234R", n2.get_reverse().name());
  {
    const auto u = n2[0];
    EXPECT_EQ((uint64_t)1234, u.id());
    EXPECT_EQ('F', u.ori());
    EXPECT_EQ(u.id(), n2.unitig_id(0));
  }

  {
    const auto u = n2[1];
    EXPECT_EQ(super_read_name::invalid_id, u.id());
    EXPECT_EQ(u.id(), n2.unitig_id(1));
  }

  static const int nb = 10;
  std::string sr;
  for(int i = 0; i < nb; ++i)
    sr += std::to_string(2 * i) + (i % 3 == 1 ? 'F' : 'R') + (i < nb - 1 ? "_" : "");
  super_read_name n3(sr);
  EXPECT_EQ((size_t)nb, n3.nb_unitigs());
  for(int i = 0; i < nb; ++i) {
    const auto u = n3[i];
    EXPECT_EQ((uint64_t)(2 * i), u.id());
    EXPECT_EQ(i % 3 == 1 ? 'F' : 'R', u.ori());
    EXPECT_EQ(u.id(), n3.unitig_id(i));
  }

  {
    const auto u = n3[nb];
    EXPECT_EQ(super_read_name::invalid_id, u.id());
    EXPECT_EQ(u.id(), n3.unitig_id(nb));
  }
}

TEST(SuperReadName, Overlap) {
  super_read_name empty("");
  super_read_name sr1("1F_2R_3F_4R");
  super_read_name sr2("4R_5F_6R");
  super_read_name sr2r("4F_5R_6F");
  super_read_name sr3("1F_2R_7F_1F_2R");
  super_read_name sr4("2R");

  EXPECT_EQ(0, empty.overlap(empty));
  EXPECT_EQ(0, empty.overlap(sr1));
  EXPECT_EQ(0, sr1.overlap(empty));
  EXPECT_EQ(1, sr1.overlap(sr2));
  EXPECT_EQ(0, sr2.overlap(sr1));
  EXPECT_EQ(0, sr1.overlap(sr2r));
  EXPECT_EQ(0, sr2r.overlap(sr1));
  EXPECT_EQ(2, sr3.overlap(sr1));
  EXPECT_EQ(1, sr3.overlap(sr4));
  EXPECT_EQ(0, sr3.overlap(sr2));
} // SuperReadName.Overlap

TEST(SuperReadName, NameOutput) {
  super_read_name    empty("");
  EXPECT_EQ("", empty.name());
  {
    std::ostringstream os; os << empty;
    EXPECT_EQ("", os.str());
  }

  super_read_name onef("54312F");
  super_read_name oner("652R");
  EXPECT_EQ("54312F", onef.name());
  EXPECT_EQ("652R", oner.name());
  {
    std::ostringstream os; os << onef;
    EXPECT_EQ("54312F", os.str());
  }
  {
    std::ostringstream os; os << oner;
    EXPECT_EQ("652R", os.str());
  }

  super_read_name many("1R_2F_65340R_123F");
  EXPECT_EQ("1R_2F_65340R_123F", many.name());
  {
    std::ostringstream os; os << many;
    EXPECT_EQ("1R_2F_65340R_123F", os.str());
  }

}

TEST(SuperReadName, Append) {
  super_read_name sr;
  EXPECT_EQ("", sr.name());

  {
    super_read_name sra;
    sr.append(sra);
    EXPECT_EQ("", sr.name());
  }
  {
    super_read_name sra("1R_2F");
    sr.append(sra);
    EXPECT_EQ("1R_2F", sr.name());
  }
  {
    super_read_name sra("1F_2R");
    sr.append(sra, 1);
    EXPECT_EQ("1R_2F_2R", sr.name());
  }
  {
    super_read_name sra("3F_4R");
    sr.append(sra, 2);
    EXPECT_EQ("1R_2F_2R", sr.name());
  }
}

TEST(SuperReadName, Prepend) {
  super_read_name sr(10);
  EXPECT_EQ("0F_0F_0F_0F_0F_0F_0F_0F_0F_0F", sr.name());
  size_t offset;

  {
    super_read_name sra("3F_4R");
    offset = sr.prepend(sra);
    EXPECT_EQ("0F_0F_0F_0F_0F_0F_0F_0F_3F_4R", sr.name());
    EXPECT_EQ((size_t)8, offset);
  }
  super_read_name sra("1F_2F_3F_4F_5F_6F_7F_8F_9F");
  offset = sr.prepend(sra, 0, offset);
  EXPECT_EQ("0F_0F_0F_0F_0F_0F_0F_0F_3F_4R", sr.name());
  EXPECT_EQ((size_t)8, offset);

  offset = sr.prepend(sra, 3, offset);
  EXPECT_EQ("0F_0F_1F_2F_3F_4F_5F_6F_3F_4R", sr.name());
  EXPECT_EQ((size_t)2, offset);

  offset = sr.prepend(sra, 7, offset);
  EXPECT_EQ("1F_2F_1F_2F_3F_4F_5F_6F_3F_4R", sr.name());
  EXPECT_EQ((size_t)0, offset);
}

void expect_sr_sequence(const std::string& res, const super_read_name& sr, const std::vector<std::string>& us, int k_len) {
  std::ostringstream os;
  sr.print_sequence(os, us, k_len);
  EXPECT_EQ(res, os.str());
}

TEST(SuperReadName, Sequence) {
  const std::vector<std::string> us = { "CGAACCTCAAGGGTTCGAT", "AAGAGTCTGTCAAGG",
                                        "GGAGTGCTCGGAAGCTGT", "GAGATCCAGCGGCTC",
                                        "GAACGACTTTAGAGTTAC" };
  const int k_len = 5;

  super_read_name sr;
  expect_sr_sequence("", sr, us, k_len);
  sr = "1F";
  expect_sr_sequence(us[1], sr, us, k_len);
  sr = "1F_3F";
  expect_sr_sequence(us[1] + us[3].substr(k_len), sr, us, k_len);

  sr = "2R";
  expect_sr_sequence(misc::rev_comp(us[2]), sr, us, k_len);
  sr = "2R_4F";
  expect_sr_sequence(misc::rev_comp(us[2]) + us[4].substr(k_len), sr, us, k_len);
  sr = "2R_4F_0R";
  expect_sr_sequence(misc::rev_comp(us[2]) + us[4].substr(k_len) + misc::rev_comp(us[0]).substr(k_len), sr, us, k_len);
} // SuperReadName.Sequence

} // empty namespace

#include <gtest/gtest.h>
#include <src_jf_aligner/super_read_name.hpp>

namespace {
TEST(SuperReadName, Parse) {
  super_read_name n1("");
  EXPECT_EQ((size_t)0, n1.nb_unitigs());
  EXPECT_EQ("", n1.get_reverse().name());
  {
    const auto u = n1[0];
    EXPECT_EQ(std::numeric_limits<unsigned long>::max(), u.id);
    EXPECT_EQ('\0', u.ori);
    EXPECT_EQ(u.id, n1.unitig_id(0));
  }

  super_read_name n2("1234F");
  EXPECT_EQ((size_t)1, n2.nb_unitigs());
  EXPECT_EQ("1234R", n2.get_reverse().name());
  {
    const auto u = n2[0];
    EXPECT_EQ((unsigned long)1234, u.id);
    EXPECT_EQ('F', u.ori);
    EXPECT_EQ(u.id, n2.unitig_id(0));
  }

  {
    const auto u = n2[1];
    EXPECT_EQ(std::numeric_limits<unsigned long>::max(), u.id);
    EXPECT_EQ('\0', u.ori);
    EXPECT_EQ(u.id, n2.unitig_id(1));
  }

  static const int nb = 10;
  std::string sr;
  for(int i = 0; i < nb; ++i)
    sr += std::to_string(2 * i) + (i % 3 == 1 ? 'F' : 'R') + (i < nb - 1 ? "_" : "");
  super_read_name n3(sr);
  EXPECT_EQ((size_t)nb, n3.nb_unitigs());
  for(int i = 0; i < nb; ++i) {
    const auto u = n3[i];
    EXPECT_EQ((unsigned long)(2 * i), u.id);
    EXPECT_EQ(i % 3 == 1 ? 'F' : 'R', u.ori);
    EXPECT_EQ(u.id, n3.unitig_id(i));
  }

  {
    const auto u = n3[nb];
    EXPECT_EQ(std::numeric_limits<unsigned long>::max(), u.id);
    EXPECT_EQ('\0', u.ori);
    EXPECT_EQ(u.id, n3.unitig_id(nb));
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

} // empty namespace

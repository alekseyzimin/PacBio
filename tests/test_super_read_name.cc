#include <gtest/gtest.h>
#include <src_jf_aligner/super_read_name.hpp>

namespace {
TEST(SuperReadName, Parse) {
  super_read_name n1("");
  EXPECT_EQ((size_t)0, n1.nb_unitigs());
  EXPECT_EQ("", n1.reverse());
  {
    const auto u = n1[0];
    EXPECT_EQ(std::numeric_limits<unsigned long>::max(), u.id);
    EXPECT_EQ('\0', u.ori);
    EXPECT_EQ(u.id, n1.unitig_id(0));
  }

  super_read_name n2("1234F");
  EXPECT_EQ((size_t)1, n2.nb_unitigs());
  EXPECT_EQ("1234R", n2.reverse());
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

} // empty namespace

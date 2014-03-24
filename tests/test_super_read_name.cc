#include <gtest/gtest.h>
#include <src_jf_aligner/super_read_name.hpp>

namespace {
TEST(SuperReadName, Parse) {
  super_read_name n1("");
  EXPECT_EQ((size_t)0, n1.nb_unitigs());
  EXPECT_EQ("", n1.reverse());
  {
    auto u = n1[0];
    EXPECT_EQ("", u.name);
    EXPECT_EQ('\0', u.ori);
  }

  super_read_name n2("1234F");
  EXPECT_EQ((size_t)1, n2.nb_unitigs());
  EXPECT_EQ("1234R", n2.reverse());
  {
    auto u = n2[0];
    EXPECT_EQ("1234", u.name);
    EXPECT_EQ('F', u.ori);
  }
  {
    auto u = n2[1];
    EXPECT_EQ("", u.name);
    EXPECT_EQ('\0', u.ori);
  }

  static const int nb = 10;
  std::string sr;
  for(int i = 0; i < nb; ++i)
    sr += std::to_string(2 * i) + (i % 3 == 1 ? 'F' : 'R') + (i < nb - 1 ? "_" : "");
  super_read_name n3(sr);
  EXPECT_EQ((size_t)nb, n3.nb_unitigs());
  for(int i = 0; i < nb; ++i) {
    auto u = n3[i];
    EXPECT_EQ(std::to_string(2 * i), u.name);
    EXPECT_EQ(i % 3 == 1 ? 'F' : 'R', u.ori);
  }
  {
    auto u = n3[nb];
    EXPECT_EQ("", u.name);
    EXPECT_EQ('\0', u.ori);
  }
}

} // empty namespace

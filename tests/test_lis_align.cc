#include <cstdlib>
#include <gtest/gtest.h>
#include <src_lis/lis_align.hpp>
#include <vector>

namespace {
struct pair {
  int first, second;
};
typedef std::vector<unsigned int> res_type;

TEST(SameSign, Int) {
  EXPECT_TRUE(lis_align::same_sign(1, 1));
  EXPECT_TRUE(lis_align::same_sign(-1, -1));
  EXPECT_FALSE(lis_align::same_sign(-1, 1));
  EXPECT_FALSE(lis_align::same_sign(1, -1));
} // SameSign.Int

bool increasing(const res_type& r, const pair* v) {
  auto it = r.cbegin();
  int  x  = v[*it].second;
  for(++it; it != r.cend(); ++it) {
    if(x >= v[*it].second)
      return false;
    x = v[*it].second;
  }
  return true;
}

TEST(LisAlign, Simple) {
  static const pair v[5] = { {1,1}, {3, 3}, {4,2}, {5,5}, {7,7} };
  auto res = lis_align::indices(v, v + 5, 5, 1);

  EXPECT_TRUE(increasing(res, v));
  EXPECT_EQ((size_t)4, res.size());
  EXPECT_EQ((unsigned int)0, res[0]);
  EXPECT_EQ((unsigned int)1, res[1]);
  EXPECT_EQ((unsigned int)3, res[2]);
  EXPECT_EQ((unsigned int)4, res[3]);
}

TEST(LisAlign, Signs) {
  {
    static const pair v[5] = { {1,-3}, {2,-2}, {3,-1}, {4,4}, {5,5} };
    auto res = lis_align::indices(v, v+ 5, 5, 1);

    EXPECT_TRUE(increasing(res, v));
    EXPECT_EQ((size_t)3, res.size());
    EXPECT_EQ((unsigned int)0, res[0]);
    EXPECT_EQ((unsigned int)1, res[1]);
    EXPECT_EQ((unsigned int)2, res[2]);
  }
  {
    static const pair v[5] = { {1,-3}, {2,-2}, {3,3}, {4,4}, {5,5} };
    auto res = lis_align::indices(v, v+ 5, 5, 1);

    EXPECT_TRUE(increasing(res, v));
    EXPECT_EQ((size_t)3, res.size());
    EXPECT_EQ((unsigned int)2, res[0]);
    EXPECT_EQ((unsigned int)3, res[1]);
    EXPECT_EQ((unsigned int)4, res[2]);
  }
}

TEST(LisAlign, Stretch) {
  static const pair v[5] = { {1, 1}, {2, 2}, {3, 3}, {10, 4} };
  auto res = lis_align::indices(v, v + 5, 5, 1);

  EXPECT_TRUE(increasing(res, v));
  EXPECT_EQ((size_t)3, res.size());
}
} // empty namespace

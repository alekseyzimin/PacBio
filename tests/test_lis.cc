#include <gtest/gtest.h>
#include <cstdlib>
#include <vector>
#include <lis.hpp>

namespace {
TEST(LIS, Increasing) {
  static const int size = 1000;
  std::vector<int> v(size);

  for(int i = 0; i < size; ++i)
    v[i] = random();

  size_t lis_size = lis::length(v.cbegin(), v.cend());
  std::vector<int> lis_seq = lis::sequence(v.cbegin(), v.cend());

  EXPECT_EQ(lis_size, lis_seq.size());
  {
    auto it = lis_seq.cbegin();
    int prev = *it;
    for(++it; it != lis_seq.cend(); prev = *it, ++it)
      EXPECT_LE(prev, *it);
  }

  std::vector<size_t> lis_ind = lis::indices(v.cbegin(), v.cend());
  EXPECT_EQ(lis_size, lis_ind.size());
  for(size_t i = 0; i < lis_size; ++i) {
    EXPECT_EQ(lis_seq[i], v[lis_ind[i]]);
    if(i > 0)
      EXPECT_LT(lis_ind[i-1], lis_ind[i]);
  }

  std::sort(v.begin(), v.end());
  std::vector<int> sorted_seq = lis::sequence(v.cbegin(), v.cend());
  EXPECT_EQ(v.size(), sorted_seq.size());
  EXPECT_EQ(v, sorted_seq);

  std::vector<int> rev_seq = lis::sequence(v.rbegin(), v.rend());
  EXPECT_EQ((size_t)1, rev_seq.size());
}
}

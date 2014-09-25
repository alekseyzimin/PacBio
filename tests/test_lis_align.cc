#include <cstdlib>
#include <gtest/gtest.h>
#include <src_lis/lis_align.hpp>
#include <vector>
#include <utility>

namespace {
// struct pair {
//   int first, second;
// };
typedef std::pair<int, int> pair;
typedef std::vector<unsigned int> res_type;

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

TEST(LisAlign, SumBuffer) {
  std::vector<int>                   numbers;
  lis_align::sum_buffer<long>        buffer(5);
  std::default_random_engine         rng;
  std::uniform_int_distribution<int> uni(-1000,1000);

  for(size_t i = 0; i < 1024; ++i)
    numbers.push_back(uni(rng));

  EXPECT_EQ((size_t)5, buffer.size());
  EXPECT_FALSE(buffer.filled());
  EXPECT_FALSE(buffer.will_be_filled());
  EXPECT_EQ(0, buffer.sum());
  for(size_t i = 0; i < numbers.size(); ++i) {
    const int test_sum = buffer.test_sum(numbers[i]);
    EXPECT_EQ(i + 1 >= buffer.size(), buffer.will_be_filled());
    buffer.push_back(numbers[i]);
    EXPECT_EQ(i + 1 >= buffer.size(), buffer.filled());
    int s = 0;
    for(size_t j = i >= buffer.size() ? i - buffer.size() + 1 : 0; j <= i; ++j)
      s += numbers[j];
    EXPECT_EQ(s, buffer.sum());
    EXPECT_EQ(test_sum, buffer.sum());
  }
}

TEST(LisAlign, Simple) {
  static const pair v[5] = { {1,1}, {3, 3}, {4,2}, {5,5}, {7,7} };
  lis_align::affine_capped accept_mer(5, 1, 1e9);
  lis_align::linear        accept_sequence(5);
  auto res = lis_align::indices(v, v + 5, 1, accept_mer, accept_sequence);

  EXPECT_TRUE(increasing(res, v));
  EXPECT_EQ((size_t)4, res.size());
  EXPECT_EQ((unsigned int)0, res[0]);
  EXPECT_EQ((unsigned int)1, res[1]);
  EXPECT_EQ((unsigned int)3, res[2]);
  EXPECT_EQ((unsigned int)4, res[3]);
}

TEST(LisAlign, StretchSimple) {
  static const pair v[4] = { {1, 1}, {2, 2}, {3, 3}, {10, 4} };
  lis_align::affine_capped accept_mer(5, 1, 1e9);
  lis_align::linear        accept_sequence(5);
  auto res = lis_align::indices(v, v + 4, 1, accept_mer, accept_sequence);

  EXPECT_TRUE(increasing(res, v));
  EXPECT_EQ((size_t)3, res.size());
}

TEST(LisAlign, StretchFullLength) {
  static const pair v[5] = { {2, 2}, {3, 3}, {13, 4}, {14, 5}, {24, 6} };
  lis_align::affine_capped accept_mer(5, 1, 1e9);
  lis_align::linear        accept_sequence(5);

  // All window stretch are satisfied but not the full length
  // stretch for 1st and 3rd case.
  {
    auto res = lis_align::indices(v, v + 3, 2, accept_mer, accept_sequence);
    EXPECT_TRUE(increasing(res, v));
    EXPECT_EQ((size_t)2, res.size());
  }

  {
    auto res = lis_align::indices(v, v + 4, 2, accept_mer, accept_sequence);
    EXPECT_TRUE(increasing(res, v));
    EXPECT_EQ((size_t)4, res.size());
  }

  {
    auto res = lis_align::indices(v, v + 5, 2, accept_mer, accept_sequence);
    EXPECT_TRUE(increasing(res, v));
    EXPECT_EQ((size_t)4, res.size());
  }
}

TEST(LisAlign, StretchLocalWindow) {
  static const pair v[5] = { {1, 1}, {2, 2}, {3, 3}, {14, 4}, {16, 5} };
  lis_align::affine_capped accept_mer(5, 1, 1e9);
  lis_align::linear        accept_sequence(5);

  {
    // Window stretch fails between {3,3} and {14,4} although full
    // length stretch would be satisfied between {1,1} and {14,4}, or
    // between {1,1} and {16,5}.
    auto res = lis_align::indices(v, v + 5, 1, accept_mer, accept_sequence);
    EXPECT_TRUE(increasing(res, v));
    EXPECT_EQ((size_t)3, res.size());
  }
}
} // empty namespace

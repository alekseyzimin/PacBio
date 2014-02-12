#include <gtest/gtest.h>
#include <lf_forward_list.hpp>

namespace {
TEST(LfForwardList, PushFront) {
  static const int size = 10000;
  lf_forward_list<int> list;
  typedef std::iterator_traits<lf_forward_list<int>::iterator>::difference_type distance_type;

  EXPECT_TRUE(list.cbegin() == list.cend());
  EXPECT_FALSE(list.cbegin() != list.cend());
  EXPECT_EQ((distance_type)0, std::distance(list.cbegin(), list.cend()));

  #pragma omp parallel for
  for(int i = 0; i < size; ++i)
    list.push_front(i);

  EXPECT_FALSE(list.cbegin() == list.cend());
  EXPECT_TRUE(list.cbegin() != list.cend());
  EXPECT_EQ((distance_type)size, std::distance(list.cbegin(), list.cend()));

  std::set<int> unique;
  for(auto it = list.cbegin(); it != list.cend(); ++it) {
    auto res = unique.insert(*it);
    EXPECT_TRUE(res.second); // Means it was not present in the set
    EXPECT_TRUE(*it >= 0 && *it < size);
  }
  EXPECT_EQ((std::set<int>::size_type)size, unique.size());

  list.clear();
  EXPECT_EQ((distance_type)0, std::distance(list.cbegin(), list.cend()));
}
}
